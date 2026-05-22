# Submodule Pinning

`superassp` vendors several upstream C/C++ libraries via git submodules. For
reproducible builds, every submodule below must be pinned to a specific commit
SHA rather than tracking a branch head.

To inspect the currently-checked-out commits, run:

```sh
git submodule status
```

| Submodule | Repository | Pinning policy |
|-----------|-----------|----------------|
| `src/SPTK` | https://github.com/sp-nitech/SPTK | Pin to a release tag (e.g. `v4.x`). |
| `src/ESTK` | https://github.com/festvox/speech_tools | Pin to a commit on `master`. |
| `src/tcl-snack` | https://github.com/scottypitcher/tcl-snack | Pin to a stable commit. |
| `src/opensmile` | https://github.com/audeering/opensmile | Pin to `v3.x`. |
| `src/Yin-Pitch-Tracking` | https://github.com/ashokfernandez/Yin-Pitch-Tracking | Pin to `master`. |
| `src/pyin` | https://github.com/aguai/pyin | Pin to `master`. |
| `src/tandem` | https://github.com/mrhuke/tandem | Pin to `master`. |
| `inst/onnx/swift-f0` | https://github.com/lars76/swift-f0 | Pin to `v0.1.1`. |

## Updating

When intentionally bumping a submodule:

1. `cd <submodule>`
2. `git fetch && git checkout <new-tag-or-sha>`
3. `cd ..` and `git add <submodule>`
4. Test thoroughly — submodule code is compiled into `superassp` and changes
   can affect numerical fidelity of `trk_*` / `lst_*` outputs.
5. Update this file with the new pinned SHA / tag.

## Cleanup

The submodule mapping in `.gitmodules` and the on-disk state may drift over
time. If `git submodule status` reports paths not present in `.gitmodules`
(e.g. `src/Parselmouth`), either remove the orphaned directory or re-add the
entry to `.gitmodules`.
