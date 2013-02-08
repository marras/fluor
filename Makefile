.SUFFIXES:

vpath %.cpp src
vpath %.cu src
vpath %.c src
vpath %.hpp include
vpath %.h include
vpath %.o obj
vpath %.oc obj

EXEC = ./fluor
EXEC_CUDA = ./fluor_CUDA
INCLUDES = -Iinclude -I/usr/include/gromacs
#CFLAGS = -g3 -Wall 
#CFLAGS = -g -pg
CFLAGS = -O2
OBJ_DIR = obj
OBJS = fluorescence.o correlatorSimple.o photonHistogram.o main.o log.o diffusion.o diffusion_obstacles.o diff_rods.o crossCorrelator.o #CUDA_random.o
OBJS_CUDA = fluorescence_c.o correlatorSimple_c.o photonHistogram_c.o main_c.o log_c.o diffusion_c.o diffusion_obstacles_c.o diff_rods_c.o crossCorrelator_c.o CUDA_random_c.o
SRC_DIR = src
#CUDA_FLAGS = -L /home/marek/NVIDIA_GPU_Computing_SDK/C/lib -l cutil_x86_64
#CUDA_FLAGS = -L ~/NVIDIA_GPU_Computing_SDK/C/lib -l cutil_i386
CUDA_FLAGS = 
#-L/usr/local/lib/gromacs -lgmx
#LIBS = -lmd -lm -lnsl -lfftw3
LIBS = -lm
#LIBS = -lm -pg	 # PG - for profiler
CC = g++
NVCC = nvcc


#Linkowanie, np.: gcc main.o -o build/hello -lstdc++
$(EXEC) : $(OBJS)
	@echo -e "---- Linking objects ----"
	cd $(OBJ_DIR); $(CC) $(OBJS) $(LIBS) $(CUDA_FLAGS) -o ../$(EXEC)
	@echo -e "\nCompilation finished successfully. You can now execute $(EXEC).\n"

#Zestaw zaleznosci + wzor kompilacji poszczegolnych czesci
fluorescence.o : fluorescence.h
diffusion.o : fluorescence.h
correlatorSimple.o : correlator.h
photonHistogram.o : photonHistogram.h
main.o: correlator.h photonHistogram.h fluorescence.h
%.o : %.cpp
	$(CC) $(INCLUDES) $(CFLAGS) -c $< -o $(OBJ_DIR)/$@
%.o : %.cu
	$(CC) $(INCLUDES) -c $< -o $(OBJ_DIR)/$@
#czyli: gcc -Iinclude -Wall -c src/plik.cpp -o obj/plik.o

##### TO SAMO W WERSJI DLA CUDY #####
$(EXEC_CUDA) : $(OBJS_CUDA)
	@echo -e "---- Linking objects ----"
	cd $(OBJ_DIR); $(NVCC) $(OBJS_CUDA) $(LIBS) $(CUDA_FLAGS) -o ../$(EXEC_CUDA) -DENABLE_GPU
	@echo -e "\nCompilation finished successfully. You can now execute $(EXEC_CUDA).\n"

#Zestaw zaleznosci + wzor kompilacji poszczegolnych czesci
fluorescence_c.o : fluorescence.h
diffusion_c.o : fluorescence.h
correlatorSimple_c.o : correlator.h
photonHistogram_c.o : photonHistogram.h
main_c.o: correlator.h photonHistogram.h fluorescence.h
%_c.o : %.cpp
	$(CC) $(INCLUDES) $(CFLAGS) -c $< -o $(OBJ_DIR)/$@ -DENABLE_GPU
%_c.o : %.cu
	$(NVCC) $(INCLUDES) -c $< -o $(OBJ_DIR)/$@ -DENABLE_GPU
#czyli: gcc -Iinclude -Wall -c src/plik.cpp -o obj/plik.o , #define ENABLE_GPU


#Delete all object + executable (in order to compile once more from scratch)
.PHONY : clean rebuild all
clean : 
	@-rm obj/*.o
	@-rm $(EXEC)
	@-rm $(EXEC_CUDA)
	@echo "*All Clean*"

rebuild : clean
rebuild : $(EXEC)

# need to write "make all" to build CUDA version too
all : $(EXEC) $(EXEC_CUDA)

cuda : $(EXEC_CUDA)

# ; - oznacza, ze nastepna komenda ma byc w tym samym shellu, \ - ze to wszystko to jest jedna linia
# @ sprawia, ze nie wypluwa na ekran "echo Compiling..." tylko samo "Compiling..." 
# $@ - target; $< - dependency
