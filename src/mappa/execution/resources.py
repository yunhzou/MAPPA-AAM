"""Process-tree memory accounting without a machine-wide scan on Linux."""
from pathlib import Path
import sys


def tree_rss(pid):
    import psutil
    if sys.platform != 'linux':
        root = psutil.Process(pid)
        return sum(p.memory_info().rss for p in [root, *root.children(recursive=True)])
    pending, visited, total = [pid], set(), 0
    while pending:
        current = pending.pop()
        if current in visited:
            continue
        visited.add(current)
        try:
            total += psutil.Process(current).memory_info().rss
            # Linux records children per thread. Include children created by
            # library threads as well as those created by the main thread.
            for thread in (Path('/proc')/str(current)/'task').iterdir():
                pending.extend(map(int, (thread/'children').read_text().split()))
        except (FileNotFoundError, ProcessLookupError, psutil.NoSuchProcess):
            pass  # the process/thread exited during this sample
    return total
