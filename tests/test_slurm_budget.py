"""Resource admission and checkpointed execution stay independent of chemistry."""
import json
from pathlib import Path
import subprocess
import sys

import pytest

from mappa.retrosynthesis.slurm_budget import SlurmBudget


@pytest.mark.parametrize('limit', [2, 3, 8, 32, 64, 100])
def test_array_admission_includes_coordinator_and_fits_small_budgets(limit):
    for workers in range(1, limit):
        budget = SlurmBudget(limit, workers)
        assert budget.concurrency * workers + 1 <= limit
        assert (budget.concurrency + 1) * workers + 1 > limit
        assert f'--array=0-999%{budget.concurrency}' in budget.arguments('0-999')


def test_budget_cannot_silently_overallocate():
    with pytest.raises(ValueError):
        SlurmBudget(4, 4)




def test_smt_allocations_are_budgeted_and_used():
    budget=SlurmBudget(3,worker_cpus=1,coordinator_cpus=1,allocation_quantum=2)
    assert budget.worker_cpus == 2
    assert budget.concurrency == 1
    assert '--cpus-per-task=2' in budget.arguments('0-63')
    with pytest.raises(ValueError):
        SlurmBudget(2,worker_cpus=1,coordinator_cpus=1,allocation_quantum=2)
