for sql in $(ls ./queries/*.sql); do
    echo "Running $sql" >&2
    python3 generate_report.py "$1" "$sql"
done