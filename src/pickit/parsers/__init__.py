"""File-format parsers for pickit: Arpeggio and IChem interaction exports.

Split out of ``io_mixin.py`` during Fase 3 (see
plan-accion-modularizacion-pickit.md), where they had been left as nested
closures inside ``analyze_files`` because separating them required turning
implicit closure state (``protein``/``ligand``/``subunit``) into explicit
parameters — a real refactor, deferred from Fase 2 on purpose to keep that
step to a same-behavior extraction.
"""
