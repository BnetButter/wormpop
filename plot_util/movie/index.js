
document.addEventListener('DOMContentLoaded', (event) => {
    // Load the SVG file
    fetch('./init_graph.svg')
        .then(response => response.text())
        .then(svg => {
            document.getElementById('svg-container').innerHTML = svg;
            const edges = [
                "egg_larva",
                "larva_adult",
                "larva_dauer",
                "dauer_larva",
                "adult_parlad",
                "egg_adult",
                "parlad_dauer",
                "larva_larvaCull",
                "larva_larvaStarve",
                "dauer_dauerCull",
                "dauer_dauerStarve",
                "parlad_parladStarve"

            ]
            const svgElement = document.querySelector('#svg-container svg');

            for (const id of edges) {
                console.log(id)
                let path = svgElement.getElementById(id).querySelector("path")
                path.setAttribute("stroke", "white")
            }

        });
});

// Function to get the playback speed
function getPlaybackSpeed() {
    const playbackSpeedInput = document.getElementById('playback-speed');
    return parseFloat(playbackSpeedInput.value);
}

function manipulateSVG(stageCount, transitionCount) {
    data = stageCount
    const svgElement = document.querySelector('#svg-container svg');
    const playbackSpeed = getPlaybackSpeed()

    const getNode = (id) => svgElement.getElementById(id).querySelector("ellipse")
    const getEdge = (id) => svgElement.getElementById(id).querySelector("polygon")

    const adultNode = getNode("adult")
    const larvaNode = getNode("larva")
    const eggNode = getNode("egg")
    const parladNode = getNode("parlad")
    const dauerNode = getNode("dauer")


    const egg_larva_edge = getEdge("egg_larva")
    const larva_adult_edge = getEdge("larva_adult")

    const larva_dauer_edge = getEdge("larva_dauer")
    const dauer_larva_edge = getEdge("dauer_larva")

    const adult_parlad_edge = getEdge("adult_parlad")

    const adult_egg_edge = getEdge("egg_adult")
    const parlad_dauer_edge = getEdge("parlad_dauer")





    let timestep = 0

    const setSize = (node, name, i) => {
        let size = Math.sqrt(parseInt(data[name][i]))
        node.setAttribute('rx', size / 3)
        node.setAttribute('ry', size / 3)
    }

    const setEdgeSize = (node, name, i) => {
        let size = Math.sqrt(parseInt(transitionCount[name][i]))
        node.setAttribute("stroke-width", size)
    }

    let execute_frame = () => {
        setSize(adultNode, "adult", timestep)
        setSize(larvaNode, "larva", timestep)
        setSize(eggNode, "egg", timestep)
        setSize(parladNode, "parlad", timestep)
        setSize(dauerNode, "dauer", timestep)

        setEdgeSize(egg_larva_edge, "egg_to_larva", timestep)
        setEdgeSize(larva_adult_edge, "larva_to_adult", timestep)
        setEdgeSize(larva_dauer_edge, "larva_to_dauer", timestep)
        setEdgeSize(dauer_larva_edge, "dauer_to_larva", timestep)
        setEdgeSize(adult_parlad_edge, "adult_to_bag", timestep)


        timestep++
    }

    setInterval(execute_frame, playbackSpeed * 100)



    // Add more manipulation code here
}


document.addEventListener('DOMContentLoaded', () => {
    document.getElementById('playback-button').addEventListener('click', processTSVFile);
});

function processTSVFile() {
    const fileInput = document.getElementById('file-upload');
    const file = fileInput.files[0];

    const stageTransition = document.getElementById("transition-upload");
    const transitionFile = stageTransition.files[0];

    if (file && transitionFile) {
        // Read the first file
        const reader1 = new FileReader();
        reader1.onload = function (e) {
            const contents1 = e.target.result;
            const stageCount = parseTSV(contents1);

            // Read the second file
            const reader2 = new FileReader();
            reader2.onload = function (e) {
                const contents2 = e.target.result;
                const transitionCount = parseTransitionTSV(contents2);
                console.log(transitionCount)
                // Call manipulateSVG with both data sets
                manipulateSVG(stageCount, transitionCount);
            };
            reader2.readAsText(transitionFile);
        };
        reader1.readAsText(file);
    } else {
        alert('Please upload both TSV files.');
    }
}


function parseTransitionTSV(content2) {
    // Split the TSV content into lines
    const lines = content2.trim().split('\n');

    // Split the first line to get headers
    const headers = lines[0].split('\t');

    // Define the headers we're interested in
    const desiredHeaders = [
        'egg_to_larva', 'egg_to_larva_mass', 'larva_to_adult', 'larva_to_adult_mass',
        'larva_to_dauer', 'larva_to_dauer_mass', 'adult_to_bag', 'adult_to_bag_mass',
        'dauer_to_larva', 'darva_to_larva_mass'
    ];

    // Create an object to store our data
    let transitions = {};

    // Iterate over each line of the TSV content, skipping the header line
    for (let i = 1; i < lines.length; i++) {
        // Split the line into its components
        const lineData = lines[i].split('\t');

        // Iterate over the desired headers and extract the corresponding data
        desiredHeaders.forEach((header, index) => {
            // Find the index of the header in the TSV file
            const headerIndex = headers.indexOf(header);

            if (headerIndex !== -1 && lineData[headerIndex] !== undefined) {
                // If the header is found and data is available, add it to the transitions object
                if (!transitions[header]) {
                    transitions[header] = [];
                }
                transitions[header].push(lineData[headerIndex]);
            }
        });
    }

    return transitions;
}

function parseTSV(tsv) {
    const lines = tsv.trim().split('\n');
    const headers = lines[0].split('\t');

    const dataIndexes = {
        larva: headers.indexOf('Number Larvae'),
        adult: headers.indexOf('Number Adults'),
        parlad: headers.indexOf('Number Parlads'),
        dauer: headers.indexOf('Number Dauer'),
        egg: headers.indexOf('Number Eggs')
    };

    const data = { larva: [], adult: [], parlad: [], dauer: [], egg: [] };

    for (let i = 1; i < lines.length; i++) {
        const row = lines[i].split('\t');
        for (const key in dataIndexes) {
            if (dataIndexes[key] !== -1) {
                data[key].push(row[dataIndexes[key]]);
            }
        }
    }

    return data;
}