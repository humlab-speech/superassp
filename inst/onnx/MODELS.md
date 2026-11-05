# ONNX models and runtime

`superassp` does not bundle neural-network weights. The ONNX model files and the
ONNX Runtime native library are downloaded on first use and cached in the R user
directory (`tools::R_user_dir("superassp", "cache")`); subsequent calls read the
cached copy with no network access. This keeps the package tarball small and
avoids shipping large binaries.

## Model licences

All weights and the inference runtime are permissively licensed (MIT), so they
may be downloaded, cached, and redistributed by downstream users and packages
without restriction beyond retaining the copyright notice.

| Model | Upstream | Licence | Hugging Face repo (pinned) | Size (approx.) |
|-------|----------|---------|-----------------------------|----------------|
| SwiftF0 | [lars76/swift-f0](https://github.com/lars76/swift-f0) (Nieradzik et al. 2025) | MIT | `FredrikKarlssonSpeech/swift-f0-onnx` @ `ab37191652a37a99cef118bc37715e20048e0fb4`, file `onnx/model.onnx` | ~400 KB |
| CREPE `tiny` | [maxrmorrison/torchcrepe](https://github.com/maxrmorrison/torchcrepe) (Kim et al. 2018) | MIT | `FredrikKarlssonSpeech/torchcrepe-onnx` @ `5768ab3ad99977bc83b15dd69a91b4d6e20987e8`, file `onnx/tiny.onnx` | ~2 MB |
| CREPE `full` | [maxrmorrison/torchcrepe](https://github.com/maxrmorrison/torchcrepe) (Kim et al. 2018) | MIT | `FredrikKarlssonSpeech/torchcrepe-onnx` @ `5768ab3ad99977bc83b15dd69a91b4d6e20987e8`, file `onnx/full.onnx` | ~85 MB |

Revisions are pinned to a commit SHA in the call site
(`R/ssff_cpp_swiftf0.R`, `R/ssff_cpp_crepe.R`) so a model cannot change under the
package without a version bump. The attribution requirement of the MIT licence is
met by this document together with the `\insertCite{}` entries in the man pages.

## ONNX Runtime

The inference runtime is ONNX Runtime v1.24.3 (MIT, Microsoft), downloaded from
GitHub releases (`https://github.com/microsoft/onnxruntime/releases`) as a
platform-specific native library (`libonnxruntime.dylib` / `libonnxruntime.so` /
`onnxruntime.dll`) and cached under
`tools::R_user_dir("superassp", "cache")/onnxruntime/`. See `R/onnxruntime.R`.

## Browser / WebAssembly

The download-and-cache mechanism above is native-R only: it relies on the
`huggingfaceR` package and on `tools::R_user_dir()` filesystem caching, and the
inference backend is a native `libonnxruntime` binary. None of these run in a
browser (WebAssembly). A browser build would instead need ONNX Runtime Web
(`onnxruntime-web`), model download via the Hugging Face `…/resolve/{revision}/{file}`
endpoints (CORS-enabled), and caching through the Cache API / IndexedDB / OPFS.
That is out of scope for this package.
