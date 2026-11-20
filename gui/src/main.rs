mod commands;

pub fn run() {
    tauri::Builder::default()
        .manage(commands::JobRegistry::default())
        .invoke_handler(tauri::generate_handler![
            commands::default_engines,
            commands::run_case,
            commands::run_case_text,
            commands::run_case_texts,
            commands::validate_case,
            commands::run_transient_text,
            commands::start_case,
            commands::poll_job,
            commands::stop_job,
            commands::read_text_file,
            commands::write_text_file,
        ])
        .run(tauri::generate_context!())
        .expect("failed to run ChuanFlow GUI");
}

fn main() {
    run();
}
