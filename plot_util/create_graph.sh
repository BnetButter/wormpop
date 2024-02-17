# Define the source DOT file and the target SVG file
SOURCE="movie/init_graph.dot"
TARGET="movie/init_graph.svg"

# Check if the source file exists
if [ -f "$SOURCE" ]; then
    # Use Graphviz's dot command to convert DOT to SVG
    dot -Tsvg "$SOURCE" -o "$TARGET"
    echo "Conversion completed: $SOURCE -> $TARGET"
else
    echo "Source file $SOURCE not found."
fi