# Historical campaign source

Retired campaign drivers and replaced visualization builders are removed from
the current checkout. Their exact source remains in [Git history](https://github.com/yunhzou/MAPPA-AAM/tree/34f6ada83adf73745d426857e088335efcd5d77a/bench/archive).
Use that revision together with the original campaign's recorded input and
environment pins for historical reproduction. Current library defaults are not
a substitute for historical settings.

The [removal manifest](../../docs/validation/repository_cleanup_20261028.json)
records each removed path, SHA256 and size, verified against that revision before
removal. The earlier [relocation manifest](../../docs/publication_cleanup_20260916.json)
retains the original report-to-archive path mapping. Measurements, witnesses and
publication artifacts remain in their established locations.

Use the [benchmark guide](../README.md) for maintained entry points.
