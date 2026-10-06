"""Score a bundled real TS guess using only the public MAPPA API."""
import json
from pathlib import Path
import numpy as np
from mappa import (
    AAMProblem, AAMSearchConfig, AtomBijection, MolecularEndpoint,
    ReactionContext, ResolvedMechanism, TransitionStateTarget, VibrationalModes,
    analyze_transition_state, ts_record,
)


def replay(path):
    data = json.loads(Path(path).read_text())
    def endpoint(key):
        item = data[key]
        return MolecularEndpoint(tuple(item['elements']),
                                 np.array(item['coordinates']),
                                 np.array(item['wbo']), label=key)
    reactant, product, guess = (endpoint(key) for key in
                                ('reactant', 'product', 'target'))
    mapping = AtomBijection.from_mapping(
        {int(a): b for a, b in data['mapping_RP'].items()},
        degree=reactant.atom_count)
    mechanism = ResolvedMechanism(
        mapping, tuple(map(tuple, data['broken'])),
        tuple(map(tuple, data['formed'])), tuple(data['core']))
    reaction = ReactionContext(
        AAMProblem(reactant, product, name=data['case']),
        AAMSearchConfig(), (mechanism,))
    target = TransitionStateTarget(
        guess, VibrationalModes(np.array(data['frequencies']),
                                np.array(data['displacements'])),
        kind='initial_guess')
    return analyze_transition_state(reaction, target)


if __name__ == '__main__':
    import signal
    def watchdog(*_):
        raise TimeoutError('TS example exceeded 300 seconds')
    signal.signal(signal.SIGALRM, watchdog)
    signal.alarm(300)
    root = Path(__file__).resolve().parents[2]
    result = replay(root / 'manuscript/animations/ts_mode_selection/inputs/recorded-case.json')
    print(json.dumps(ts_record(result), indent=2))
