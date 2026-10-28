from rebenchmark.relocate_index import relocate


def test_relocation_changes_only_explicit_mount_prefix():
    original = dict(path='/old/data/a', sha256='abc', count=4,
                    nested=[('/old/data', '/old/database/a', '/other/a')])
    changed = relocate(original, '/old/data', '/new/data')
    assert changed == dict(path='/new/data/a', sha256='abc', count=4,
                          nested=[('/new/data', '/old/database/a', '/other/a')])
    assert original['path'] == '/old/data/a'


def test_saved_member_identity_is_unchanged():
    member = dict(archive='/old/a.tar', name='finalized/00000.pkl.gz',
                  offset=512, bytes=1234, sha256='input-hash')
    moved = relocate(member, '/old', '/new')
    assert moved.pop('archive') == '/new/a.tar'
    assert moved == {k: v for k, v in member.items() if k != 'archive'}
