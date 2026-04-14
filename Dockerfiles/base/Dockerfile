# ~~~~~~~~~~~~~~~~~~~~SLIM BASE ICEFEM IMAGE DOCKERFILE~~~~~~~~~~~~~~~~~~~~
FROM --platform=linux/amd64 ubuntu:latest AS builder
LABEL maintainer="Sarah Wells-Moran swellsmo@uchicago.edu"

# CHANGED: cleaned up/combined the installs
RUN apt-get -y update && \
    apt-get -y install apt-utils build-essential gfortran libopenmpi-dev wget unzip \
        firefox gedit gcc g++ gdb cmake-curses-gui mesa-common-dev mesa-utils freeglut3-dev \
        ninja-build pkgconf libc6-amd64-cross libstdc++6 cmake git curl xxd xdg-utils libnss3 && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

ENV rcFile=/root/.bashrc LIBRARY_DIR=/usr/lib MKL_INSTALL_FOLDER=/usr/lib/Intel

RUN touch /root/.bashrc && \
    mkdir ${MKL_INSTALL_FOLDER} && \
    echo "export LIBRARY_DIR=${LIBRARY_DIR}" >> "$rcFile"

# ── Intel oneAPI Base Toolkit ──────────────────────────────────────
WORKDIR ${MKL_INSTALL_FOLDER}
ADD https://registrationcenter-download.intel.com/akdlm/IRC_NAS/992857b9-624c-45de-9701-f6445d845359/l_BaseKit_p_2023.2.0.49397.sh ${MKL_INSTALL_FOLDER}

# CHANGED: Combined all these into one layer - it takes up less space.
RUN /bin/bash ./l_BaseKit_p_2023.2.0.49397.sh -a -s --eula accept --install-dir $MKL_INSTALL_FOLDER && \
    rm l_BaseKit_p_2023.2.0.49397.sh && \
    chmod +x ${MKL_INSTALL_FOLDER}/setvars.sh && \
    echo ". ${MKL_INSTALL_FOLDER}/setvars.sh intel64" >> "$rcFile"

# CHANGED: FYI env vars set by sourcing a script don't persist across RUN layers — only Docker ENV directives do. This is probably what caused the "No such file or directory" for petscrules/petscconf.h.
ENV MKLROOT=/usr/lib/Intel/mkl/latest
ENV I_MPI_ONEAPI_ROOT=/usr/lib/Intel/mpi/latest
ENV I_MPI_ROOT=/usr/lib/Intel/mpi/latest
ENV LD_LIBRARY_PATH=${MKLROOT}/lib/intel64:${I_MPI_ONEAPI_ROOT}/lib:${LD_LIBRARY_PATH}
ENV PATH=${I_MPI_ONEAPI_ROOT}/bin:${MKLROOT}/bin:${PATH}
ENV CMAKE_PREFIX_PATH=${MKLROOT}:${I_MPI_ONEAPI_ROOT}:${CMAKE_PREFIX_PATH}

# ── PETSc (debug build) ───────────────────────────────────────────
WORKDIR /
ENV PETSC_DIR=/usr/lib/petsc
# CHANGED: Added the git clone, pinning to v3.22.3
RUN git clone -b v3.22.3 --depth 1 https://gitlab.com/petsc/petsc.git ${PETSC_DIR}
WORKDIR ${PETSC_DIR}

ENV PETSC_ARCH=debug
RUN ./configure PETSC_ARCH=debug \
        --with-scalapack-include=$MKLROOT/include \
        --with-scalapack-lib="-L$MKLROOT/lib/intel64 -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64" \
        --with-mpi-dir=$I_MPI_ONEAPI_ROOT \
        --with-mkl=$MKLROOT \
        --with-blaslapack-dir=$MKLROOT \
        --with-mkl_pardiso-dir=$MKLROOT \
        --with-mkl_cpardiso-dir=$MKLROOT \
        --download-hdf5 \
        --download-mumps \
        --with-clanguage=cxx \
        --with-debugging=1 && \
    make PETSC_DIR=${PETSC_DIR} PETSC_ARCH=debug all && \
    make PETSC_DIR=${PETSC_DIR} PETSC_ARCH=debug check

# ── PETSc (optimised build) ───────────────────────────────────────
ENV PETSC_ARCH=optimised
RUN ./configure PETSC_ARCH=optimised \
        --with-scalapack-include=$MKLROOT/include \
        --with-scalapack-lib="-L$MKLROOT/lib/intel64 -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64" \
        --with-mpi-dir=$I_MPI_ONEAPI_ROOT \
        --with-mkl=$MKLROOT \
        --with-blaslapack-dir=$MKLROOT \
        --with-mkl_pardiso-dir=$MKLROOT \
        --with-mkl_cpardiso-dir=$MKLROOT \
        --download-hdf5 \
        --download-mumps \
        --with-clanguage=cxx \
        --with-debugging=0 \
        --COPTFLAGS='-O3' --FOPTFLAGS='-O3' --CXXOPTFLAGS='-O3' && \
    make PETSC_DIR=${PETSC_DIR} PETSC_ARCH=optimised all && \
    make PETSC_DIR=${PETSC_DIR} PETSC_ARCH=optimised check

# ── Third-party libraries ─────────────────────────────────────────
WORKDIR /

ENV RAPIDJSON_DIR=${LIBRARY_DIR}/rapidjson
ENV EIGEN3_Dir=${LIBRARY_DIR}/eigen
ENV HIGHFIVE_DIR=${LIBRARY_DIR}/HighFive
# CHANGED: VTK_DIR must point to the build directory (where VTKConfig.cmake gets generated), not the source directory. CMake couldn't find VTK otherwise.
ENV VTK_DIR=${LIBRARY_DIR}/VTK/build

# CHANGED: FYI Docker's ADD with a URL doesn't clone repos.
RUN git clone --depth 1 https://github.com/Tencent/rapidjson.git  ${RAPIDJSON_DIR} && \
    git clone --depth 1 https://gitlab.com/libeigen/eigen.git      ${EIGEN3_Dir}    && \
    git clone --depth 1 https://github.com/BlueBrain/HighFive.git  ${HIGHFIVE_DIR}  && \
    git clone --depth 1 --branch v9.3.1 https://gitlab.kitware.com/vtk/vtk.git ${LIBRARY_DIR}/VTK && \
    wget -O ${LIBRARY_DIR}/igafem.zip https://sourceforge.net/projects/cmcodes/files/latest/download && \
    mkdir ${LIBRARY_DIR}/igafem && \
    unzip ${LIBRARY_DIR}/igafem.zip -d ${LIBRARY_DIR}/igafem && \
    rm ${LIBRARY_DIR}/igafem.zip

WORKDIR ${LIBRARY_DIR}/VTK
RUN mkdir -p build && cd build && \
    cmake -GNinja ../ && \
    cmake --build .

# ── Ice_FEM application ───────────────────────────────────────────
WORKDIR /
RUN git clone --depth 1 https://github.com/swellsmo/Ice_FEM_2025a.git /Ice_FEM_2025a
WORKDIR /Ice_FEM_2025a
RUN cmake . && \
    make -j$(nproc) && \
    mkdir -p Results

# ── Clean up build artifacts before the runtime
RUN rm -rf ${PETSC_DIR}/.git \
           ${PETSC_DIR}/src \
           ${PETSC_DIR}/debug/externalpackages \
           ${PETSC_DIR}/optimised/externalpackages \
           ${RAPIDJSON_DIR}/.git \
           ${EIGEN3_Dir}/.git \
           ${HIGHFIVE_DIR}/.git \
           ${LIBRARY_DIR}/VTK/.git \
           ${LIBRARY_DIR}/VTK/Common ${LIBRARY_DIR}/VTK/Filters ${LIBRARY_DIR}/VTK/IO \
           ${LIBRARY_DIR}/VTK/Rendering ${LIBRARY_DIR}/VTK/CMake \
           /Ice_FEM_2025a/.git \
           /Ice_FEM_2025a/Documentation \
           ${LIBRARY_DIR}/igafem/C_files_win \
           ${LIBRARY_DIR}/igafem/compile_win.m ${LIBRARY_DIR}/igafem/compile.m

FROM --platform=linux/amd64 ubuntu:latest AS runtime

RUN apt-get -y update && \
    apt-get -y install --no-install-recommends \
        build-essential gfortran cmake ninja-build \
        libopenmpi-dev libstdc++6 libc6 \
        libgomp1 libopengl0 libglx0 ca-certificates wget && \
    apt-get clean && \
    apt-get autoremove && \
    rm -rf /var/lib/apt/lists/*

# ── Copy Intel MKL + MPI runtime libs only ────────────────────────
# We only need the MKL shared libraries and Intel MPI runtime, and third-party libraries
COPY --from=builder /usr/lib/Intel/mkl/latest/lib     /usr/lib/Intel/mkl/latest/lib
COPY --from=builder /usr/lib/Intel/mkl/latest/include /usr/lib/Intel/mkl/latest/include
COPY --from=builder /usr/lib/Intel/mpi/latest/lib     /usr/lib/Intel/mpi/latest/lib
COPY --from=builder /usr/lib/Intel/mpi/latest/bin     /usr/lib/Intel/mpi/latest/bin
COPY --from=builder /usr/lib/Intel/mpi/latest/include /usr/lib/Intel/mpi/latest/include
COPY --from=builder /usr/lib/petsc     /usr/lib/petsc
COPY --from=builder /usr/lib/rapidjson /usr/lib/rapidjson
COPY --from=builder /usr/lib/eigen     /usr/lib/eigen
COPY --from=builder /usr/lib/HighFive  /usr/lib/HighFive
COPY --from=builder /usr/lib/VTK/build /usr/lib/VTK/build
COPY --from=builder /usr/lib/igafem    /usr/lib/igafem

# ── Copy Ice_FEM application ─────────────────────────────────────
COPY --from=builder /Ice_FEM_2025a /Ice_FEM_2025a

# ── Environment variables ─────────────────────────────────────────
ENV LIBRARY_DIR=/usr/lib
ENV MKL_INSTALL_FOLDER=/usr/lib/Intel
ENV MKLROOT=/usr/lib/Intel/mkl/latest
ENV I_MPI_ONEAPI_ROOT=/usr/lib/Intel/mpi/latest
ENV I_MPI_ROOT=/usr/lib/Intel/mpi/latest
ENV PETSC_DIR=/usr/lib/petsc
ENV PETSC_ARCH=optimised
ENV RAPIDJSON_DIR=/usr/lib/rapidjson
ENV EIGEN3_Dir=/usr/lib/eigen
ENV HIGHFIVE_DIR=/usr/lib/HighFive
ENV VTK_DIR=/usr/lib/VTK/build
ENV LD_LIBRARY_PATH=${MKLROOT}/lib/intel64:${I_MPI_ONEAPI_ROOT}/lib
ENV PATH=${I_MPI_ONEAPI_ROOT}/bin:${MKLROOT}/bin:${PATH}
ENV CMAKE_PREFIX_PATH=${MKLROOT}:${I_MPI_ONEAPI_ROOT}
ENV MATLAB_DIR=/home/matlab

RUN touch /root/.bashrc && \
    echo "export PETSC_DIR=${PETSC_DIR}" >> /root/.bashrc && \
    echo "export PETSC_ARCH=${PETSC_ARCH}" >> /root/.bashrc && \
    echo "export VTK_DIR=${VTK_DIR}" >> /root/.bashrc && \
    echo "export LIBRARY_DIR=${LIBRARY_DIR}" >> /root/.bashrc && \
    echo "PS1='🐳 \[\e[1;32m\]\u@ice-fem\[\e[m\]:\[\e[1;34m\]\w\[\e[m\]\\$ '" >> /root/.bashrc

WORKDIR /Ice_FEM_2025a
CMD ["/bin/bash"]
