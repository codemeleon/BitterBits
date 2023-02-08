julia -e "using Pkg; \
	Pkg.add(\"DataFrames\"); \
	Pkg.add(\"Turing\"); \
	Pkg.add(\"Statistics\"); \
	Pkg.add(\"Flux\"); \
	Pkg.add(\"Distributions\"); \
	Pkg.add(\"CSV\"); \
	Pkg.add(\"ArviZ\"); \
	Pkg.add(\"Glob\"); \
	Pkg.add(\"FilePaths\"); \
	Pkg.add(\"Gadfly\"); \
	Pkg.add(\"KernelDensity\"); \
	Pkg.add(\"StatsPlots\"); \
	Pkg.add(\"Random\"); \
	Pkg.add(\"RDatasets\"); \
	Pkg.add(\"Clustering\"); \
	Pkg.add(\"Cleaner\"); \
	Pkg.add(\"TableTransforms\"); \
	Pkg.add(\"Bio\"); \
	Pkg.add(\"PlutoExtras\"); \
	Pkg.add(\"IJulia\"); \
	"

# Pkg.add(\"SciML\"); \ Need to explore later
