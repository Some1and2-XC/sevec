# Sevec
### A fast, pure Rust implementation of the segmented array data structure.

## Purpose:
The purpose of the [`Sevec`] data structure is to allow for fast array splitting and copy operations.
If an application needs to handle a large amount of data.

This library is written to use as few allocations as possible, the style of code in this library is a lot more similar to C code than traditional rust code.

## Unsafe:
This library uses a good deal of unsafe code. This is for performance.
If concerns about the soundness of the unsafe exist, the following command can be used to validate:

```sh
MIRIFLAGS=-Zmiri-strict-provenance rustup run nightly cargo miri test
```
