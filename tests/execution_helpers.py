"""Small importable operations for execution integration tests."""
import os
import time


def value(payload, workspace):
    allocation = bytearray(payload.get('allocate', 0))
    time.sleep(payload.get('sleep', 0))
    if payload.get('fail'):
        raise ValueError('deliberate test failure')
    if payload.get('exit'):
        os._exit(7)
    return payload.get('value', 1)


def total(payload, workspace):
    return sum(workspace.dependency(key) for key in payload)
