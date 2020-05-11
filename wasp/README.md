# WASP

This is the Web ASembly Port for BABY-GROOT.

> this is work in progress...

To check out the current web version of the app that is produced by this repo, go to [willrowe.net/baby-groot](https://willrowe.net/baby-groot)

Usage instructions are available on the baby-groot app.

## Running locally

To browserify the js, build the WASM binary and run the development server:

``` bash
make dev
```

Then navigate to: [http://localhost:3434/](http://localhost:3434/)

## Issues / TODO

* can't load large index
  * runs out of memory when allocating
* premature terminations aren't graceful
  * they just print to the console, no notifications for user
  * in some cases, the application doesn't shut down
