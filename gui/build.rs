fn main() {
    // No beforeBuildCommand: the frontend is pre-built and committed under
    // frontend/dist, so a missing tauri CLI is not a problem.
    tauri_build::build()
}
