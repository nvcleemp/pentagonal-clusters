
SOURCES = pentagon_partition.c filter_valid_clusters.c count_all_clusters.c\
          has_six_cluster.c appearances_of_clusters.c min_edge_count.c\
          Makefile COPYRIGHT LICENSE README.md

all: build/pentagon_partition build/filter_valid_clusters build/count_all_clusters\
     build/has_six_cluster build/appearances_of_clusters build/min_edge_count

clean:
	rm -rf build
	rm -rf dist

build/pentagon_partition: pentagon_partition.c
	mkdir -p build
	cc -o $@ -O4 -Wall $^

build/filter_valid_clusters: filter_valid_clusters.c
	mkdir -p build
	cc -o $@ -O4 -Wall $^

build/count_all_clusters: count_all_clusters.c
	mkdir -p build
	cc -o $@ -O4 -Wall $^

build/appearances_of_clusters: appearances_of_clusters.c
	mkdir -p build
	cc -o $@ -O4 -Wall $^

build/has_six_cluster: has_six_cluster.c
	mkdir -p build
	cc -o $@ -O4 -Wall $^

build/min_edge_count: min_edge_count.c
	mkdir -p build
	cc -o $@ -O4 -Wall $^

sources: dist/pentagonal-cluster-sources.zip dist/pentagonal-cluster-sources.tar.gz

dist/pentagonal-cluster-sources.zip: $(SOURCES)
	mkdir -p dist
	zip dist/pentagonal-cluster-sources $(SOURCES)

dist/pentagonal-cluster-sources.tar.gz: $(SOURCES)
	mkdir -p dist
	tar czf dist/pentagonal-cluster-sources.tar.gz $(SOURCES)
