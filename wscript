APPNAME = 'test-project'
VERSION = '1.0.0'

srcdir = '.'
blddir = 'build'
def set_options(opt):
    opt.tool_options('compiler_cxx')
    pass

def configure(conf):
    #conf.load('g++-4.9')
    conf.check_tool('compiler_cxx')
    conf.check_cxx(lib = 'blas')
    conf.check_cxx(lib = 'lapack')
    conf.check_cxx(lib = 'gfortran')
    conf.check_cxx(lib = 'omp')
    conf.check_cxx(lib = 'math')
    conf.check_cxx(lib = 'libc++')
    pass

def build(bld):
    bld(features = 'cxx cprogram',
        source = 'src/transition_matrix.cpp src/data.cpp src/ctm_data_refresh.cpp src/wright_fisher_estimater.cpp src/q_function.cpp src/WF_main.cpp src/boundary.cpp',
        include = ['include','Matrix/include'],
        lib = ['lapack','blas','gfortran'],
	#cflags = ['-fopenmp'],
        cxxflag = ['-std=c++14'],
	cxxflags = ['-Ofast','-fopenmp','-mavx','-march=native'],
        #cxxflags = ['-Og','-fopenmp'],
	linkflag = ['-Ofast','-fopenmp','-mavx','-march=native'],
        #linkflag = ['-Og','-fopenmp'],
        target = 'WF_main')
    pass

def shutdown(bld):
    pass
    
