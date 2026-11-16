"""Shared state storage retains ordinary immutable sequence semantics."""
import copy
import pickle
import pytest
from mappa.mappa_core.core.search.pairs import SharedPairs


def test_partition_independent_equality_hashing_and_indexing():
    values=((1,4),(2,5),(3,6))
    left=SharedPairs((values[:1],values[1:]),3)
    right=SharedPairs((values,),3)
    assert left==right==values
    assert values==left
    assert hash(left)==hash(right)==hash(values)
    assert {left:'saved'}[right]=='saved'
    assert dict(left)==dict(values)
    assert left[-1]==values[-1]
    assert left[1:]==values[1:]
    with pytest.raises(IndexError):
        left[3]
    assert copy.deepcopy(left) is left
    assert pickle.loads(pickle.dumps(left))==left
