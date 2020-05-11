/* global FileReader */

const { Readable } = require('readable-stream')
const toBuffer = require('typedarray-to-buffer')

class FileReadStream extends Readable {
  constructor (file, opts = {}) {
    super(opts)

    // save the read offset
    this._offset = 0
    this._ready = false
    this._file = file
    this._size = file.size
    this._chunkSize = opts.chunkSize || Math.max(this._size / 1000, 200 * 1024)

    // create the reader
    const reader = new FileReader()

    reader.onload = () => {
      // get the data chunk
      this.push(toBuffer(reader.result))
    }
    reader.onerror = () => {
      this.emit('error', reader.error)
    }

    this.reader = reader

    // generate the header blocks that we will send as part of the initial payload
    this._generateHeaderBlocks(file, opts, (err, blocks) => {
      // if we encountered an error, emit it
      if (err) {
        return this.emit('error', err)
      }

      // push the header blocks out to the stream
      if (Array.isArray(blocks)) {
        blocks.forEach(block => this.push(block))
      }

      this._ready = true
      this.emit('_ready')
    })
  }

  _generateHeaderBlocks (file, opts, callback) {
    callback(null, [])
  }

  _read () {
    if (!this._ready) {
      this.once('_ready', this._read.bind(this))
      return
    }

    const startOffset = this._offset
    let endOffset = this._offset + this._chunkSize
    if (endOffset > this._size) endOffset = this._size

    if (startOffset === this._size) {
      this.destroy()
      this.push(null)
      return
    }

    this.reader.readAsArrayBuffer(this._file.slice(startOffset, endOffset))

    // update the stream offset
    this._offset = endOffset
  }

  destroy () {
    this._file = null
    if (this.reader) {
      this.reader.onload = null
      this.reader.onerror = null
      try { this.reader.abort() } catch (e) {};
    }
    this.reader = null
  }
}

module.exports = FileReadStream
