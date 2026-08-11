"""
ZTR chunk compression codec

ZTR chunk data is not compressed once but through a *chain* of filters, each
layer prefixing its own format byte. A typical SMP4 chunk reads::

    zlib -> RLE -> FOLLOW1 -> 16TO8 -> DELTA2 -> raw

Decoding therefore loops until the leading byte is ``RAW``. Every filter stores
its own parameters in the stream (the RLE guard byte, the delta level), so an
encoded chain can be replayed exactly on write, reproducing the original bytes.

Ported from Staden io_lib's ``compression.c``, which is what SnapGene uses. Its
zlib layer is ``deflateInit2(level=1, windowBits=15, memLevel=8,
Z_HUFFMAN_ONLY)``; matching those settings is what makes round-trips
byte-identical.
"""

import logging
import struct
import zlib
from typing import Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)

MAGIC = b"\xaeZTR\r\n\x1a\n"

RAW = 0
RLE = 1
ZLIB = 2
DELTA1 = 64
DELTA2 = 65
DELTA4 = 66
SHRINK_16TO8 = 70
SHRINK_32TO8 = 71
FOLLOW1 = 72

FORMAT_NAMES = {
    RAW: "RAW", RLE: "RLE", ZLIB: "ZLIB", 3: "XRLE", 4: "XRLE2",
    DELTA1: "DELTA1", DELTA2: "DELTA2", DELTA4: "DELTA4",
    67: "DDELTA1", 68: "DDELTA2", 69: "DDELTA4",
    SHRINK_16TO8: "16TO8", SHRINK_32TO8: "32TO8", FOLLOW1: "FOLLOW1",
    73: "CHEB445", 74: "ICHEB", 75: "LOG2", 77: "STHUFF",
    79: "QSHIFT", 80: "TSHIFT",
}

# io_lib compresses with Z_HUFFMAN_ONLY at level 1; anything else produces a
# different byte stream for identical input.
_ZLIB_LEVEL = 1
_ZLIB_STRATEGY = zlib.Z_HUFFMAN_ONLY

_MAX_LAYERS = 16  # guards against a malformed stream looping forever


class ZtrFormatError(Exception):
    """A ZTR filter this codec does not implement."""


def _signed(byte: int) -> int:
    return byte - 256 if byte > 127 else byte


# =============================================================================
# DECODERS
# =============================================================================


def _un_zlib(data: bytes) -> bytes:
    return zlib.decompress(data[5:])


def _un_rle(data: bytes) -> bytes:
    out_len = struct.unpack("<I", data[1:5])[0]
    guard = data[5]
    out = bytearray()
    i = 6
    n = len(data)
    while len(out) < out_len and i < n:
        byte = data[i]
        i += 1
        if byte != guard:
            out.append(byte)
            continue
        count = data[i]
        i += 1
        if count == 0:
            # 'guard 0' escapes a literal guard byte in the input
            out.append(guard)
        else:
            out.extend(bytes([data[i]]) * count)
            i += 1
    return bytes(out)


def _un_delta1(data: bytes) -> bytes:
    level = data[1]
    body = data[2:]
    out = bytearray(len(body))
    u1 = u2 = u3 = 0
    for i, byte in enumerate(body):
        if level == 1:
            z = u1
        elif level == 2:
            z = 2 * u1 - u2
            u2 = u1
        else:
            z = 3 * u1 - 3 * u2 + u3
            u3 = u2
            u2 = u1
        u1 = _signed((_signed(byte) + z) & 0xFF)
        out[i] = u1 & 0xFF
    return bytes(out)


def _un_delta2(data: bytes) -> bytes:
    level = data[1]
    body = data[2:]
    out = bytearray(len(body))
    u1 = u2 = u3 = 0
    for i in range(0, len(body) - 1, 2):
        if level == 1:
            z = u1
        elif level == 2:
            z = 2 * u1 - u2
            u2 = u1
        else:
            z = 3 * u1 - 3 * u2 + u3
            u3 = u2
            u2 = u1
        u1 = (((body[i] << 8) | body[i + 1]) + z) & 0xFFFFFFFF
        out[i] = (u1 >> 8) & 0xFF
        out[i + 1] = u1 & 0xFF
    return bytes(out)


def _un_delta4(data: bytes) -> bytes:
    level = data[1]
    body = data[4:]  # 2 dummy bytes pad the header to a 4-byte boundary
    out = bytearray(len(body))
    u1 = u2 = u3 = 0
    for i in range(0, len(body) - 3, 4):
        if level == 1:
            z = u1
        elif level == 2:
            z = 2 * u1 - u2
            u2 = u1
        else:
            z = 3 * u1 - 3 * u2 + u3
            u3 = u2
            u2 = u1
        value = (body[i] << 24) | (body[i + 1] << 16) | (body[i + 2] << 8) | body[i + 3]
        u1 = (value + z) & 0xFFFFFFFF
        out[i] = (u1 >> 24) & 0xFF
        out[i + 1] = (u1 >> 16) & 0xFF
        out[i + 2] = (u1 >> 8) & 0xFF
        out[i + 3] = u1 & 0xFF
    return bytes(out)


def _un_16to8(data: bytes) -> bytes:
    out = bytearray()
    i = 1
    n = len(data)
    while i < n:
        value = _signed(data[i])
        if value >= 0:
            out += b"\x00" + data[i : i + 1]
            i += 1
        elif value != -128:
            out += b"\xff" + data[i : i + 1]
            i += 1
        else:
            # -128 escapes a value that did not fit in 8 bits
            i += 1
            out += data[i : i + 2]
            i += 2
    return bytes(out)


def _un_32to8(data: bytes) -> bytes:
    out = bytearray()
    i = 1
    n = len(data)
    while i < n:
        value = _signed(data[i])
        if value != -128:
            fill = b"\xff" if value < 0 else b"\x00"
            out += fill * 3 + data[i : i + 1]
            i += 1
        else:
            i += 1
            out += data[i : i + 4]
            i += 4
    return bytes(out)


def _un_follow1(data: bytes) -> bytes:
    successor = data[1:257]
    body = data[257:]
    if not body:
        return b""
    out = bytearray([body[0]])
    for i in range(1, len(body)):
        out.append((successor[out[i - 1]] - _signed(body[i])) & 0xFF)
    return bytes(out)


_DECODERS = {
    RLE: _un_rle,
    ZLIB: _un_zlib,
    DELTA1: _un_delta1,
    DELTA2: _un_delta2,
    DELTA4: _un_delta4,
    SHRINK_16TO8: _un_16to8,
    SHRINK_32TO8: _un_32to8,
    FOLLOW1: _un_follow1,
}


# =============================================================================
# ENCODERS
# =============================================================================


def _do_zlib(data: bytes) -> bytes:
    compressor = zlib.compressobj(_ZLIB_LEVEL, zlib.DEFLATED, 15, 8, _ZLIB_STRATEGY)
    body = compressor.compress(data) + compressor.flush()
    return bytes([ZLIB]) + struct.pack("<I", len(data)) + body


def _pick_guard(data: bytes) -> int:
    """Least frequent byte, lowest value winning ties — as io_lib chooses it."""
    counts = [0] * 256
    for byte in data:
        counts[byte] += 1
    return min(range(256), key=lambda i: (counts[i], i))


def _do_rle(data: bytes, guard: Optional[int] = None) -> bytes:
    if guard is None:
        guard = _pick_guard(data)
    out = bytearray()
    i = 0
    n = len(data)
    while i < n:
        k = i
        while k < n and data[i] == data[k]:
            if k - i == 255:
                break
            k += 1
        if k - i < 4:
            # Runs of 1-3 are cheaper stored literally
            while True:
                if data[i] == guard:
                    out += bytes([guard, 0])
                else:
                    out.append(data[i])
                i += 1
                if k < i + 1:
                    break
        else:
            out += bytes([guard, k - i, data[i]])
            i = k
    return bytes([RLE]) + struct.pack("<I", n) + bytes([guard]) + bytes(out)


def _do_delta1(data: bytes, level: int) -> bytes:
    out = bytearray(len(data))
    u1 = u2 = u3 = 0
    for i in range(len(data)):
        if level == 1:
            z = u1
        elif level == 2:
            z = 2 * u1 - u2
            u2 = u1
        else:
            z = 3 * u1 - 3 * u2 + u3
            u3 = u2
            u2 = u1
        u1 = _signed(data[i])
        out[i] = (u1 - z) & 0xFF
    return bytes([DELTA1, level]) + bytes(out)


def _do_delta2(data: bytes, level: int) -> bytes:
    out = bytearray(len(data))
    u1 = u2 = u3 = 0
    for i in range(0, len(data) - 1, 2):
        if level == 1:
            z = u1
        elif level == 2:
            z = 2 * u1 - u2
            u2 = u1
        else:
            z = 3 * u1 - 3 * u2 + u3
            u3 = u2
            u2 = u1
        u1 = (data[i] << 8) | data[i + 1]
        delta = (u1 - z) & 0xFFFF
        out[i] = (delta >> 8) & 0xFF
        out[i + 1] = delta & 0xFF
    return bytes([DELTA2, level]) + bytes(out)


def _do_delta4(data: bytes, level: int) -> bytes:
    out = bytearray(len(data))
    u1 = u2 = u3 = 0
    for i in range(0, len(data) - 3, 4):
        if level == 1:
            z = u1
        elif level == 2:
            z = 2 * u1 - u2
            u2 = u1
        else:
            z = 3 * u1 - 3 * u2 + u3
            u3 = u2
            u2 = u1
        u1 = (data[i] << 24) | (data[i + 1] << 16) | (data[i + 2] << 8) | data[i + 3]
        delta = (u1 - z) & 0xFFFFFFFF
        out[i] = (delta >> 24) & 0xFF
        out[i + 1] = (delta >> 16) & 0xFF
        out[i + 2] = (delta >> 8) & 0xFF
        out[i + 3] = delta & 0xFF
    return bytes([DELTA4, level, 0, 0]) + bytes(out)


def _do_16to8(data: bytes) -> bytes:
    out = bytearray([SHRINK_16TO8])
    for i in range(0, len(data) - 1, 2):
        value = struct.unpack(">h", data[i : i + 2])[0]
        if -127 <= value <= 127:
            out.append(value & 0xFF)
        else:
            out.append(0x80)
            out += data[i : i + 2]
    return bytes(out)


def _do_32to8(data: bytes) -> bytes:
    out = bytearray([SHRINK_32TO8])
    for i in range(0, len(data) - 3, 4):
        value = struct.unpack(">i", data[i : i + 4])[0]
        if -127 <= value <= 127:
            out.append(value & 0xFF)
        else:
            out.append(0x80)
            out += data[i : i + 4]
    return bytes(out)


def _do_follow1(data: bytes) -> bytes:
    successor = bytearray(256)
    best = [0] * 256
    seen: Dict[int, int] = {}
    for i in range(len(data) - 1):
        current, following = data[i], data[i + 1]
        key = (current << 8) | following
        count = seen.get(key, 0) + 1
        seen[key] = count
        if count > best[current]:
            best[current] = count
            successor[current] = following
    out = bytearray([FOLLOW1])
    out += successor
    if data:
        out.append(data[0])
        for i in range(1, len(data)):
            out.append((successor[data[i - 1]] - _signed(data[i])) & 0xFF)
    return bytes(out)


# =============================================================================
# CHAIN
# =============================================================================


def unpack(data: bytes) -> Tuple[bytes, List[Dict]]:
    """Decode a chunk down to raw bytes, returning ``(raw, chain)``.

    ``chain`` records each layer outermost-first along with the parameters
    needed to reapply it, so ``pack()`` can reproduce the input byte for byte.
    Raises :class:`ZtrFormatError` for a filter this codec does not implement.
    """
    chain: List[Dict] = []
    for _ in range(_MAX_LAYERS):
        if not data or data[0] == RAW:
            return data, chain
        form = data[0]
        decoder = _DECODERS.get(form)
        if decoder is None:
            raise ZtrFormatError(
                f"unsupported ZTR filter {FORMAT_NAMES.get(form, form)}"
            )
        layer: Dict = {"format": form}
        if form == RLE:
            layer["guard"] = data[5]
        elif form in (DELTA1, DELTA2, DELTA4):
            layer["level"] = data[1]
        chain.append(layer)
        data = decoder(data)
    raise ZtrFormatError("ZTR filter chain too deep")


def pack(raw: bytes, chain: List[Dict]) -> bytes:
    """Re-apply a chain recorded by :func:`unpack`, innermost layer first."""
    data = raw
    for layer in reversed(chain):
        form = layer.get("format")
        if form == ZLIB:
            data = _do_zlib(data)
        elif form == RLE:
            data = _do_rle(data, layer.get("guard"))
        elif form == DELTA1:
            data = _do_delta1(data, layer.get("level", 1))
        elif form == DELTA2:
            data = _do_delta2(data, layer.get("level", 1))
        elif form == DELTA4:
            data = _do_delta4(data, layer.get("level", 1))
        elif form == SHRINK_16TO8:
            data = _do_16to8(data)
        elif form == SHRINK_32TO8:
            data = _do_32to8(data)
        elif form == FOLLOW1:
            data = _do_follow1(data)
        else:
            raise ZtrFormatError(
                f"cannot re-encode ZTR filter {FORMAT_NAMES.get(form, form)}"
            )
    return data
