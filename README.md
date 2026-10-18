# tclcluspack

Tcl wrapper around the ClusPack clustering library. Package `tclcluspack` 0.1 provides one command,
`Tclcluspack`, used by the Ordali alignment editor. It builds against Tcl 8.6 and Tcl 9 from the same sources.

```tcl
package require tclcluspack
Tclcluspack $alignment -dt alignment -cm bionj -nbc secator
```

`CHANGES.txt` (in French) lists the bugs fixed in 2026.

## Build

TEA (autoconf) build, out of the source tree:

```sh
mkdir build-tcl8 && cd build-tcl8
../configure --with-tcl=/path/to/tcl8.6/lib && make     # -> libtclcluspack0.1.so
mkdir build-tcl9 && cd build-tcl9
../configure --with-tcl=/path/to/tcl9/lib && make       # -> libtcl9tclcluspack0.1.so
make install
```

`lib/` holds prebuilt static Rmath libraries (from R) for linux-x86_64, macosx-arm64, macosx-x86_64
and win32-x86_64; `configure` picks the one for the build platform. On Windows, build with MinGW/MSYS2.
`configure` is already generated; run `autoconf` only after editing `configure.ac`.

Compile flags: `-DCLUSPACK_VERBOSE` enables debug output, `-DCLUSPACK_SYSTEM_MALLOC` uses
`malloc`/`free` instead of the Tcl allocator (for AddressSanitizer or valgrind).

## Tests

```sh
tests/run.sh build-tcl8                                  # every case of tests/cases.tcl, one process each
TCLSH=tclsh9.0 tests/run.sh build-tcl9                   # same under Tcl 9
tclsh8.6 tests/stress.tcl build-tcl8 100                 # repeated calls, reports memory growth
```

## License

GPL-3.0-or-later, Copyright (c) Luc Moulinier. See [LICENSE](LICENSE).
The clustering code comes from ClusPack; R's Rmath library is GPL-2.0-or-later.
