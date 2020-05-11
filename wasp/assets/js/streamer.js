'use strict'

/* ===========================================================================
   Libs
   =========================================================================== */
var FileReadStream = require('filestream/read')
var zlib = require('zlib')
var peek = require('peek-stream')
const { Transform } = require('stream')

/* ===========================================================================
   HTML ID tags
   =========================================================================== */
const $fastqUploader = document.getElementById('fastqUploader')
const $fastqFileName = document.getElementById('fastqFileName')
const $fastqSelecter = document.getElementById('fastqSelecter')
const $noFile = document.getElementById('noFile')
const $progressBar = document.getElementById('progress-bar')

/* ===========================================================================
   FASTQ handling
   =========================================================================== */
// FileDragHover is used to cancel event and hover styling
function FileDragHover(e) {
    e.stopPropagation()
    e.preventDefault()
    e.target.className = e.type == 'dragover' ? 'hover' : ''
}

// FileSelectHandler processes any files that are added
function FileSelectHandler(e) {
    FileDragHover(e)

    // fetch FileList object
    var files = e.target.files || e.dataTransfer.files
    if (files.length != 0) {
        // process all the FASTQ files
        var fastqList = []
        var fileNames = ''
        for (var i = 0, f;
            (f = files[i]); i++) {
            // record the filename for the user
            fileNames += f.name.replace('C:\\fakepath\\', '') + ','

            // add the file for GROOT
            fastqList.push(f)
        }

        // show the user the filenames
        $fastqSelecter.classList.add('active')
        $fastqFileName.textContent = 'selected file(s): '
        $noFile.innerHTML = fileNames

        // pass the FASTQ file list to WASM
        getFiles(fastqList)
        iconUpdate('inputIcon')
    } else {
        $fastqSelecter.classList.remove('active')
        $noFile.textContent = 'no files :('
    }
}

/* ===========================================================================
   FASTQ parsing (adapted from: https://blog.luizirber.org/static/sourmash-wasm/app.js)
   =========================================================================== */
const resetProgress = () => {
    $progressBar.style.transform = 'translateX(-100%)'
    $progressBar.style.opacity = '0'
}

function formatBytes(bytes) {
    if (bytes < 1024) return bytes + ' Bytes'
    else if (bytes < 1048576) return (bytes / 1024).toFixed(2) + ' KB'
    else if (bytes < 1073741824) return (bytes / 1048576).toFixed(2) + ' MB'
    else return (bytes / 1073741824).toFixed(2) + ' GB'
}

function isGzip(data) {
    return data[0] === 31 && data[1] === 139
}

function newParser() {
    return peek(function(data, swap) {
        if (isGzip(data)) return swap(null, new zlib.Unzip())
        else {
            return swap(
                null,
                new Transform({
                    transform(chunk, encoding, callback) {
                        this.push(chunk)
                        callback()
                    }
                })
            )
        }
    })
}

// this is the exported function - fastqStreamer - it is called by WASM when GROOT is ready to processes FASTQ data
module.exports = async function(fileArr) {
    var files = fileArr[0]

    // get the progress bar ready
    resetProgress()
    $progressBar.style.opacity = '1.0'

    // process the file list
    for (let file of files) {
        await streamFile(file)
        console.log('ended stream: ', file.name)
        resetProgress()
    }

    // closeFASTQchan is a close down signal, sent to WASM once all the FASTQ data has been sent
    console.log('closing the WASM FASTQ channel')
    closeFASTQchan()
}

// streamFile is called for each file sent to the fastqStreamer function
async function streamFile(file) {
    return new Promise(resolve => {
        // set up the data stream and check errors
        let stream = new FileReadStream(file)
        stream.reader.onerror = e => {
            console.error(e)
        }
        let fileSize = formatBytes(stream._size)
        let loadedFile = 0
        let lastUpdate = 0
        console.log('started stream: ', file.name)
        statusUpdate('status', '> processed 0 / ' + fileSize)

        // monitor progress
        stream.reader.onprogress = data => {
            if (data.loaded == data.total && loadedFile < stream._size) {
                loadedFile += data.loaded
                let percentLoaded = Math.round((loadedFile / file.size) * 100)
                $progressBar.style.transform = `translateX(${-(100 - percentLoaded)}%)`
                if (percentLoaded % 2 == 0 && percentLoaded != lastUpdate) {
                    statusUpdate(
                        'status',
                        '> processed ' + formatBytes(loadedFile) + ' / ' + fileSize
                    )
                    lastUpdate = percentLoaded
                }
            }
        }

        // set up the parser
        let parser = newParser()
        parser.on('data', function(data) {
            // munchFASTQ is linked to WASM and used to send the data to Go
            munchFASTQ(data, data.length)
        })

        // pipe the data stream into the parser
        stream.pipe(parser)
        stream.on('end', () => resolve())
        stream.on('error', error => reject(error))
    })
}

/* ===========================================================================
   GROOT set up (graph and index loading)
   =========================================================================== */
// getGraphs gets the groot graphs ready for loading
function getGraphs(graphURL) {
    var reader = new FileReader()
    fetch(graphURL)
        .then(function(response) {
            if (!response.ok) {
                statusUpdate('status', '> could not download groot graphs!')
            }
            return response.blob()
        })
        .then(data => {
            reader.readAsArrayBuffer(data)
        })
        .catch(function(error) {
            console.log(error)
        })
    reader.onload = ev => {
        var raw_data = new Uint8Array(
            ev.target.result,
            0,
            ev.target.result.byteLength
        )
        if (raw_data === null) {
            statusUpdate('status', '> could not download groot graphs!')
        } else {
            loadGraphs(graphURL, raw_data, reader.result.byteLength)
        }
    }
}

// getLSHforest gets the index ready for loading
function getLSHforest(lshfURL) {
    var reader = new FileReader()
    fetch(lshfURL)
        .then(function(response) {
            if (!response.ok) {
                statusUpdate('status', '> could not download index!')
            }
            return response.blob()
        })
        .then(data => {
            reader.readAsArrayBuffer(data)
        })
        .catch(function(error) {
            console.log(error)
        })
    reader.onload = ev => {
        var raw_data = new Uint8Array(
            ev.target.result,
            0,
            ev.target.result.byteLength
        )
        if (raw_data === null) {
            statusUpdate('status', '> could not download index!')
        } else {
            loadIndex(lshfURL, raw_data, reader.result.byteLength)
        }
    }
}

/* ===========================================================================
   Boot the app
   =========================================================================== */
const startApplication = () => {
    // setup the page
    window.onload = function() {
        window.ontouchmove = function() {
            return false
        }
        window.onorientationchange = function() {
            document.body.scrollTop = 0
        }

        // launch webassembly
        if (WebAssembly) {
            // WebAssembly.instantiateStreaming is not currently available in Safari
            if (WebAssembly && !WebAssembly.instantiateStreaming) {
                // polyfill
                WebAssembly.instantiateStreaming = async(resp, importObject) => {
                    const source = await (await resp).arrayBuffer()
                    return await WebAssembly.instantiate(source, importObject)
                }
            }
            const go = new Go()
            WebAssembly.instantiateStreaming(
                fetch('baby-groot.wasm'),
                go.importObject
            ).then(result => {
                go.run(result.instance)

                // set up the fastq listener
                $fastqUploader.addEventListener('change', FileSelectHandler, false)

                // load the graphs and index file
                //getGraphs('assets/groot-files/betalactam/groot.gg')
                //getLSHforest('assets/groot-files/betalactam/groot.lshe')

                getGraphs('assets/groot-files/index/groot.gg')
                getLSHforest('assets/groot-files/index/groot.lshe')
                    // getGraphs('assets/groot-files/dummy-db/groot.gg')
                    // getLSHforest('assets/groot-files/dummy-db/groot.lshe')
                toggleDiv('spinner')
                statusUpdate('status', '> GROOT is ready!')
            })
        } else {
            toggleDiv('spinner')
            console.log('WebAssembly is not supported in this browser')
            statusUpdate('status', '> please get a more recent browser!')
        }
    }

    // listen out for index selection -- TODO: I don't think this isn't listening to the right selecter
    // document.getElementById("indexSelecter").addEventListener('click', function() {
    //  var selectedIndex = document.getElementById("indexSelecter").selectedOptions;
    //  getIndex(selectedIndex[0].value);
    // });
}

startApplication()