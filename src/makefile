# Makefile

# �v���O�������ƃI�u�W�F�N�g�t�@�C����
program = rlr
objs    = scene.o vectormath.o modelholder.o mtseq.o

# ��`�σ}�N���̍Ē�`
CC     = x86_64-w64-mingw32-g++ #gcc
CXX    = x86_64-w64-mingw32-g++
CFLAGS = -Wall -O2 -funroll-loops
LDFLAGS = -static-libstdc++ -static-libgcc

# �T�t�B�b�N�X���[���K�p�Ώۂ̊g���q�̒�`
.SUFFIXES: .cxx .o

# �v���C�}���^�[�Q�b�g
$(program): $(objs)
	$(CC) $(LDFLAGS) -o $(program) $^

# �T�t�B�b�N�X���[��
.cxx.o:
	$(CXX) $(CFLAGS) -c $<

# �t�@�C���폜�p�^�[�Q�b�g
.PHONY: clean
clean:
	$(RM) $(program) $(objs)


scene.o       : scene.cxx scene.h bvh.h vectormath.h modelholder.h kdtree.h zmx.h
vectormath.o  : vectormath.cxx vectormath.h
modelholder.o : modelholder.cxx modelholder.h
mtseq.o       : mtseq.cxx mtseq.h