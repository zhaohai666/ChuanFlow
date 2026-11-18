// Minimal Tauri 2 IPC bridge, so the GUI needs no `npm install`.
// Tauri 2 injects window.__TAURI_INTERNALS__ with an `invoke` function;
// @tauri-apps/api/core simply forwards to it, and so do we.
export async function invoke(cmd, args = {}) {
  const internals = globalThis.__TAURI_INTERNALS__;
  if (!internals || typeof internals.invoke !== 'function') {
    throw new Error('Tauri IPC 不可用：请在桌面窗口内运行本界面');
  }
  return internals.invoke(cmd, args);
}
