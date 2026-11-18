mod commands;

pub fn run() {
    tauri::Builder::default()
        .invoke_handler(tauri::generate_handler![
            commands::run_case,
            commands::run_case_text,
            commands::validate_case,
        ])
        .run(tauri::generate_context!())
        .expect("failed to run ChuanFlow GUI");
}

fn main() {
    run();
}
