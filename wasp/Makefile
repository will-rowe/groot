# index parameters
DB='arg-annot'
WINDOW_SIZE=100
K_SIZE=51
SKETCH_SIZE=30
NUM_PART=8
MAX_K=4
MAX_SPAN=30
CT=0.99
THREADS=8

# output files
BINARY='baby-groot.wasm'
INDEX='assets/groot-files/index'

all: build index

dev: build server

build:
	npm run build --prefix assets/js/
	GOOS=js GOARCH=wasm go build -o ${BINARY} -ldflags "-s -w" ./cmd

index:
	go run ../main.go get -d ${DB}
	go run ../main.go index -m ${DB}* -i ${INDEX} -w ${WINDOW_SIZE} -k ${K_SIZE} -s ${SKETCH_SIZE} -x ${NUM_PART} -y ${MAX_K} --maxSketchSpan ${MAX_SPAN} -p ${THREADS}
	rm -r ${DB}*

server:
	go run dev-server.go

clean:
	rm ${BINARY}
	rm -r ${INDEX}
	rm groot.log