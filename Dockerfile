FROM debian:stretch-slim
# (c) Roel Kluin r.kluin@nki.nl GPL-2.0

COPY install/nkiBRCA_1.00.tar.gz /tmp

RUN BWA_VERSION=0.7.17 \
  && HTSLIB_VERSION=1.10.2 \
  && SAMTOOLS_VERSION=1.10 \
  && BEDTOOLS_VERSION=2.29.2 \
  && R_VERSION=3.5.1 \
#
# day before bioconductor 3.8 released. This depracated data.frame to IRanges, required for cghseg
  && BUILD_DATE=2018-10-30 \
#
# update certificates, apt-utils warnings are safe to ignore.
  && echo 'debconf debconf/frontend select Noninteractive' | debconf-set-selections \
  && export TERM=linux \
  && apt-get update && apt-get upgrade -y \
  && apt-get install -y --no-install-recommends --no-install-suggests ca-certificates \
  && update-ca-certificates \
#
# a minimal system for R and samtools
  && apt-get install -y --no-install-recommends \
    libblas-dev libbz2-1.0 libcurl3 libicu57 libjpeg62-turbo libopenblas-dev \
    libpangocairo-1.0-0 libpcre3 libpng16-16 libreadline7 libtiff5 liblzma5 \
    zlib1g fonts-texgyre gsfonts locales libblas3 libgomp1 \
#
  && echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen \
  && locale-gen en_US.utf8 \
  && /usr/sbin/update-locale LANG=en_US.UTF-8 \
  && export LC_ALL=en_US.UTF-8 \
  && export LANG=en_US.UTF-8 \
#
# prevent error when installing openjdk-8-jre-headless
  && mkdir -p /usr/share/man/man1 \
# tools for building
  && BUILDDEPS="default-jdk libbz2-dev libcairo2-dev libcurl4-openssl-dev\
       libpango1.0-dev libjpeg-dev libicu-dev libpcre3-dev libpng-dev\
       libreadline-dev libtiff5-dev liblzma-dev libx11-dev libxt-dev\
       perl tcl8.6-dev tk8.6-dev texinfo texlive-extra-utils gfortran\
       texlive-fonts-recommended texlive-fonts-extra texlive-latex-recommended\
       x11proto-core-dev xauth xfonts-base xvfb zlib1g-dev wget curl\
       make gcc g++ libc-dev build-essential autoconf zip unzip\
       bzip2 file dpkg-dev git" \
#
  && apt-get install -y --no-install-recommends --no-install-suggests $BUILDDEPS \
#
#
################## R #######################
#
  && cd /tmp \
#
  && wget https://cran.r-project.org/src/base/R-3/R-${R_VERSION}.tar.gz \
#
  && tar -xf R-${R_VERSION}.tar.gz \
#
  && cd /tmp/R-${R_VERSION} \
#
  && R_PAPERSIZE=A4 \
    R_BATCHSAVE="--no-save --no-restore" \
    PERL=/usr/bin/perl \
    R_UNZIPCMD=/usr/bin/unzip \
    LIBnn=lib \
    AWK=/usr/bin/awk \
    CFLAGS="-g -O2 -fstack-protector-strong -Wformat -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2 -g" \
    CXXFLAGS="-g -O2 -fstack-protector-strong -Wformat -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2 -g" \
  ./configure --with-readline --with-blas --disable-nls \
#
  && make \
#
  && mkdir -p /usr/local/lib/R/lib/ \
#
  && make install \
#
## Add a default CRAN mirror
  && echo "options(repos = c(CRAN = 'https://cran.rstudio.com/'), download.file.method = 'libcurl')" >> /usr/local/lib/R/etc/Rprofile.site \
#
## Add a library directory (for user-installed packages)
  && mkdir -p /usr/local/lib/R/site-library \
  && chown root:staff /usr/local/lib/R/site-library \
  && chmod g+wx /usr/local/lib/R/site-library \
#
## Fix library path
  && echo "R_LIBS_USER='/usr/local/lib/R/site-library'" >> /usr/local/lib/R/etc/Renviron \
  && echo "R_LIBS=\${R_LIBS-'/usr/local/lib/R/site-library:/usr/local/lib/R/library:/usr/lib/R/library'}" >> /usr/local/lib/R/etc/Renviron \
#
## Clean up from R source install
  && cd /tmp && rm -rf R-${R_VERSION}* \
#
#
################## bwa #######################
#
#
# dangling head warnings are safe to ignore.
  && git clone https://github.com/lh3/bwa.git \
#
  && cd bwa; git checkout v${BWA_VERSION} && make \
#
  && mv bwa /usr/local/bin/ \
  && cd /tmp && rm -r bwa \
#
#
################## htslib #######################
#
  && git clone https://github.com/samtools/htslib.git \
#
  && cd htslib; git checkout ${HTSLIB_VERSION} \
        && autoheader \
        && autoconf -Wno-syntax \
        && ./configure --prefix=/usr/local \
        && make \
        && make install \
#
  && cd /tmp && rm -r htslib \
#
#
################## samtools #######################
#
  && git clone https://github.com/samtools/samtools.git \
#
  && cd samtools; git checkout ${SAMTOOLS_VERSION} \
        && autoheader \
        && autoconf -Wno-syntax \
        && ./configure --prefix=/usr/local --without-curses \
        && make \
        && make install \
  && cd /tmp && rm -r samtools \
#
  && rm -rf /usr/local/share/man \
#
#
################## bedtools #######################
#
  && git clone https://github.com/arq5x/bedtools2.git \
#
  && cd bedtools2; git fetch --all --tags; git checkout v${BEDTOOLS_VERSION} \
  && make \
  && cp -f "bin/"* "/usr/local/bin" \
#
  && cd /tmp && rm -r bedtools2 \
#
################## R packages #######################
#
## install packages from date-locked MRAN snapshot of CRAN
  && [ -z "$BUILD_DATE" ] && BUILD_DATE=$(TZ="America/Los_Angeles" date -I) || true \
  && MRAN=https://mran.microsoft.com/snapshot/${BUILD_DATE} \
  && echo MRAN=$MRAN >> /etc/environment \
  && export MRAN=$MRAN \
  && echo "options(repos = c(CRAN='$MRAN'), download.file.method = 'libcurl')" >> /usr/local/lib/R/etc/Rprofile.site \
#
## install packages per one package should ensure depenendencies are installed as well
 && Rscript -e "install.packages('BiocManager', INSTALL_opts = c('--no-html', '--no-help', '--no-docs', '--no-demo', '--no-multiarch', '--no-data'));\
 BiocManager::install('BiocVersion', INSTALL_opts = c('--no-html', '--no-help', '--no-docs', '--no-demo', '--no-multiarch', '--no-data'), update=FALSE);\
 BiocManager::install('BiocGenerics', INSTALL_opts = c('--no-html', '--no-help', '--no-docs', '--no-demo', '--no-multiarch', '--no-data'), update=FALSE);\
 BiocManager::install('S4Vectors', INSTALL_opts = c('--no-html', '--no-help', '--no-docs', '--no-demo', '--no-multiarch', '--no-data'), update=FALSE);\
 BiocManager::install('BiocVersion', INSTALL_opts = c('--no-html', '--no-help', '--no-docs', '--no-demo', '--no-multiarch', '--no-data'), update=FALSE);\
 BiocManager::install('IRanges', INSTALL_opts = c('--no-html', '--no-help', '--no-docs', '--no-demo', '--no-multiarch', '--no-data'), update=FALSE);\
 install.packages('gtools', INSTALL_opts = c('--no-html', '--no-help', '--no-docs', '--no-demo', '--no-multiarch', '--no-data'));\
 install.packages('readxl', INSTALL_opts = c('--no-html', '--no-help', '--no-docs', '--no-demo', '--no-multiarch', '--no-data'));\
 install.packages('RODBC', INSTALL_opts = c('--no-html', '--no-help', '--no-docs', '--no-demo', '--no-multiarch', '--no-data')); \
 install.packages('rjson', INSTALL_opts = c('--no-html', '--no-help', '--no-docs', '--no-demo', '--no-multiarch', '--no-data')); \
 install.packages('optparse', INSTALL_opts = c('--no-html', '--no-help', '--no-docs', '--no-demo', '--no-multiarch', '--no-data')); \
 install.packages('pamr', INSTALL_opts = c('--no-html', '--no-help', '--no-docs', '--no-demo', '--no-multiarch', '--no-data')); \
 install.packages('zoo', INSTALL_opts = c('--no-html', '--no-help', '--no-docs', '--no-demo', '--no-multiarch', '--no-data'));"
#
#
################## cghseg #######################
#
RUN wget https://cran.r-project.org/src/contrib/Archive/cghseg/cghseg_1.0.5.tar.gz \
  && R CMD INSTALL --no-html --no-docs --no-help --no-demo --no-multiarch --no-data cghseg_1.0.5.tar.gz \
  && R CMD INSTALL --no-html --no-docs --no-help --no-demo --no-multiarch --no-data /tmp/nkiBRCA_1.00.tar.gz \
  && rm cghseg_* /tmp/nkiBRCA* \
#

################## nkiBRCA #######################
#
		&& wget http://ccb.nki.nl/software/nkibrca/nkiBRCA_1.00.tar.gz \
		&& R CMD INSTALL --no-html --no-docs --no-help --no-demo --no-multiarch --no-data nkiBRCA_1.00.tar.gz \
		&& rm nkiBRCA_1.00* \
#
		
################## clean up #########################
# packages with + appended are kept
  && apt-get remove --purge -y $BUILDDEPS libblas3+ libgomp1+ \
#
  && apt-get autoremove -y \
  && apt-get autoclean -y \
#
  && rm -rf /var/lib/apt \
#
################## command #########################
#
  && mkdir /app

COPY *.sh *.R ref /app/

CMD ["/app/align_count_and_classify.sh"]

