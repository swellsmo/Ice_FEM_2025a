#~~~~~~~~~~~~~~~~~~~~~~~~~~MATLAB ADDITION TO ICEFEM ~~~~~~~~~~~~~~~~~~~~~~~~~~~
FROM --platform=linux/amd64 icefem:base AS builder
LABEL maintainer="Sarah Wells-Moran swellsmo@uchicago.edu"

# No clue if I have to redefine these again but can't hurt ¯\_(ツ)_/¯  
ENV rcFile=/root/.bashrc LIBRARY_DIR=/usr/lib MKL_INSTALL_FOLDER=/usr/lib/Intel
ENV MKLROOT=/usr/lib/Intel/mkl/latest
ENV I_MPI_ONEAPI_ROOT=/usr/lib/Intel/mpi/latest
ENV I_MPI_ROOT=/usr/lib/Intel/mpi/latest
ENV LD_LIBRARY_PATH=${MKLROOT}/lib/intel64:${I_MPI_ONEAPI_ROOT}/lib:${LD_LIBRARY_PATH}
ENV PATH=${I_MPI_ONEAPI_ROOT}/bin:${MKLROOT}/bin:${PATH}
ENV CMAKE_PREFIX_PATH=${MKLROOT}:${I_MPI_ONEAPI_ROOT}:${CMAKE_PREFIX_PATH}
ENV PETSC_DIR=/usr/lib/petsc
ENV PETSC_ARCH=optimised
ENV RAPIDJSON_DIR=${LIBRARY_DIR}/rapidjson
ENV EIGEN3_Dir=${LIBRARY_DIR}/eigen
ENV HIGHFIVE_DIR=${LIBRARY_DIR}/HighFive
ENV VTK_DIR=${LIBRARY_DIR}/VTK/build

# Stage 2: RUNTIME — slim image with only what's needed to run
FROM --platform=linux/amd64 mathworks/matlab:r2024a AS runtime

# Add any necessary MATLAB toolboxes here, separated by spaces
ARG ADDITIONAL_PRODUCTS="Parallel_Computing_Toolbox Image_Processing_Toolbox"

# Default user for matlab image is matlab. Set user to root to install runtime dependencies
USER 0
RUN apt-get -y update && \
    apt-get -y install --no-install-recommends \
        build-essential gfortran cmake ninja-build \
        libopenmpi-dev libstdc++6 libc6 \
        libgomp1 libopengl0 libglx0 ca-certificates wget && \
    apt-get clean && \
    apt-get autoremove && \
    rm -rf /var/lib/apt/lists/*

# Install Matlab Package Manager to download toolboxes specified in ADDITIONAL_PRODUCTS
WORKDIR /tmp
USER matlab
RUN wget -q https://www.mathworks.com/mpm/glnxa64/mpm \
    && chmod +x mpm \
    && EXISTING_MATLAB_LOCATION=$(dirname $(dirname $(readlink -f $(which matlab)))) \
    && sudo HOME=${HOME} ./mpm install \
        --destination=${EXISTING_MATLAB_LOCATION} \
        --release=R2024a \
        --products ${ADDITIONAL_PRODUCTS} \
    || (echo "MPM Installation Failure. See below for more information:" && cat /tmp/mathworks_root.log && false) \
    && sudo rm -rf mpm /tmp/mathworks_root.log

WORKDIR /home/matlab

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

CMD ["/bin/bash"]
