"""SHA-256 hex digests of files and byte strings, as manifests and receipts record them."""
import hashlib
from pathlib import Path


def sha256(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def sha256_bytes(data):
    return hashlib.sha256(data).hexdigest()
