"""Explicit, trusted-only conversion of saved fragment/1 catalogues for decoding.

This is an experiment input conversion, never a runtime compatibility fallback.
It removes the retired self-tolerance filter while retaining saved binary
generators, discovery constraints, placements, histories and original weights.
It does not rerun matching or claim the catalogue came from the current matcher.
"""
from dataclasses import asdict, replace
import gc
import gzip
import hashlib
import io
import json

from mappa import AAMSearchConfig
from mappa.mappa_core.core.scheduler.bundles import Blob, BundleWriter
from mappa.mappa_core.core.scheduler.results import CatalogueResult
from mappa.mappa_core.core.symmetry.structure import BinarySymmetry
from mappa.mappa_core.utils.io.checkpoints import _CheckpointUnpickler


def search_config(record):
    """Remove only the explicitly retired filter; retain its value in provenance."""
    record = dict(record)
    record.pop('self_tolerance')
    return AAMSearchConfig(**record)


class _RetiredConstraint:
    """Unpickling carrier; discarded with the retired self-tolerance filter."""


class _RetiredConfig:
    pass


class _BinaryConversion:
    @classmethod
    def _from_binary(cls, atoms, elements, anchors, threshold, tolerance,
                     generators, constraints, weights, edges, projected):
        # These generators already describe the full binary graph. The old
        # constraints restricted their joint action by WBO self-tolerance.
        return BinarySymmetry._from_binary(
            atoms, elements, anchors, threshold, generators, weights, edges)


class _MigrationUnpickler(_CheckpointUnpickler):
    def find_class(self, module, name):
        if module == 'mappa.mappa_core.core.symmetry.relaxed':
            return {'RelaxedAutomorphism': _BinaryConversion,
                    'BondConstraint': _RetiredConstraint}[name]
        if (module, name) == ('mappa.mappa_core.domain', 'AAMSearchConfig'):
            return _RetiredConfig
        return super().find_class(module, name)


def convert_blob(reference, expected_config):
    """Check the immutable range before explicitly changing its symmetry policy."""
    reference = Blob(**reference) if isinstance(reference, dict) else reference
    with open(reference.path, 'rb') as stream:
        stream.seek(reference.offset)
        raw = stream.read(reference.size)
    if len(raw) != reference.size or hashlib.sha256(raw).hexdigest() != reference.sha256:
        raise ValueError('Historical catalogue checksum differs')
    with gzip.GzipFile(fileobj=io.BytesIO(raw)) as stream:
        saved = _MigrationUnpickler(stream).load()
    if (not isinstance(saved, CatalogueResult)
            or json.loads(json.dumps(vars(saved.search_config))) != expected_config):
        raise ValueError('Historical catalogue configuration differs')
    catalogue = saved.catalogue
    for family in catalogue.families:
        family.validate_representative(catalogue.problem)
    for fragment in catalogue._fragments.values():
        if not fragment.relation.allows_discovery(fragment.witness):
            raise ValueError('Converted fragment lost its discovery witness')
    return replace(saved, search_config=search_config(expected_config))


def convert_rows(rows, directory, expected_config):
    """Write aggregate bundles; original references and policy remain auditable."""
    writer = BundleWriter(directory, None)
    converted, evidence = [], []
    for name, row in rows:
        saved = convert_blob(row['catalogue'], expected_config)
        if saved.case != name:
            raise ValueError('Catalogue direction differs from its ledger')
        reference = asdict(writer.put(saved))
        catalogue = saved.catalogue
        converted.append((name, dict(row, catalogue=reference)))
        evidence.append(dict(case=name, original=row['catalogue'], converted=reference,
            fragments=len(catalogue._fragments), families=len(catalogue.families),
            branches=catalogue.branch_count, paths=catalogue.path_count))
        del saved, catalogue
        gc.collect()
        if len(converted) % 100 == 0:
            print(f'Converted {len(converted)}/{len(rows)} saved catalogues', flush=True)
    return converted, evidence
