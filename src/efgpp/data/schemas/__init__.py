"""Pandera contracts for every standardized table the Data layer writes.

User input is checked by the adapters (which report every problem with counts); these
schemas guard the *outputs* so downstream layers can rely on them.
"""
