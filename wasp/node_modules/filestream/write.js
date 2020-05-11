/* global File */

const { Writable } = require('readable-stream')
const toBuffer = require('typedarray-to-buffer')

class FileWriteStream extends Writable {
  constructor (callback, opts) {
    // inherit writable
    super(Object.assign({ decodeStrings: false }, opts))

    // when the stream finishes create a file
    this.on('finish', this._generateFile.bind(this))

    // create the internal buffers storage
    this._buffers = []
    this._bytesreceived = 0
    this.callback = callback
    this.type = (opts || {}).type
  }

  _createFile () {
    // if we have no buffers, then abort any processing
    if (this._buffers.length === 0) {
      return
    }

    return new File(this._buffers, '', {
      type: this.type || ''
    })
  }

  _generateFile () {
    const file = this._createFile()

    if (file) {
      if (typeof this.callback === 'function') {
        this.callback(file)
      }

      this.emit('file', file)
    }

    // reset the buffers and counters
    this._buffers = []
    this._bytesreceived = 0
  }

  _preprocess (data, callback) {
    // pass through the data
    callback(null, data)
  }

  _write (chunk, encoding, callback) {
    const data = Buffer.isBuffer(chunk) ? chunk : toBuffer(chunk)
    const writeStream = this

    this._preprocess(data, (err, processed) => {
      if (err) {
        return callback(err)
      }

      // if the incoming data has been passed through,
      // then add to the bytes received buffer
      if (processed) {
        writeStream._bytesreceived += processed.length
        writeStream._buffers.push(processed)
        writeStream.emit('progress', writeStream._bytesreceived)
      }

      callback()
    })
  }
}

module.exports = FileWriteStream
