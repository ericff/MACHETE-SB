FROM nathankw/centos6
#Comes with Python v2.7.10, JRE v1.8.0_91, and R v3.2.3.
#Also comes with Perl v5.10.1, which was installed as part of the "Development Tools" package.
#The directories /srv/src and /srv/software are created in the base image. 
MAINTAINER Nathaniel Watson <nathankw@stanford.edu>
RUN yum install -y bc 
#INSTALL TBB (Threading Building Blocks) from Intel.
# Can be used when intalling Bowtie1 and Bowtie2 in order to use TBB over pthreads for managing parallel processes.
#RUN mkdir /srv/src/TBB && \
#	cd /srv/src/TBB && \
#	wget https://www.threadingbuildingblocks.org/sites/default/files/software_releases/source/tbb44_20160128oss_src_0.tgz && \
#	tar -zxf tbb44_20160128oss_src_0.tgz && \
#	cd tbb44_20160128oss && \
#	gmake && \
#	. build/linux_intel64_gcc_cc4.4.7_libc2.12_kernel4.1.19_release/tbbvars.sh

#INSTALL Bowtie1.1.1
RUN mkdir /srv/src/Bowtie1 && \
	cd /srv/src/Bowtie1 && \
	wget https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.1.2/bowtie-1.1.2-src.zip && \
	unzip bowtie-1.1.2-src.zip && \
	cd bowtie-1.1.2 && \
#	make WITH_TBB=1 && \
	make && \
	make install
	
#INSTALL Bowtie2 2.2.8
RUN mkdir /srv/src/Bowtie2 && \
	cd /srv/src/Bowtie2 && \
	wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.8/bowtie2-2.2.8-source.zip && \
	unzip bowtie2-2.2.8-source.zip && \
	cd bowtie2-2.2.8 && \
#	make WITH_TBB=1 && \
	make && \
	make install
#INSTALL R packages chron and data.table.
# chron: Chronological Objects which can Handle Dates and Times (https://cran.r-project.org/web/packages/chron/index.html)
#	data.table: Extension of data.frame (https://cran.r-project.org/web/packages/data.table/data.table.pdf)
RUN wget https://cran.r-project.org/src/contrib/chron_2.3-47.tar.gz && \
		R CMD INSTALL chron_2.3-47.tar.gz && \
		wget https://cran.r-project.org/src/contrib/data.table_1.9.6.tar.gz && \	
		R CMD INSTALL data.table_1.9.6.tar.gz

#INSTALL samtools/1.3. Needed for knife.
RUN mkdir /srv/src/samtools && \
	cd /srv/src/samtools && \
	wget https://github.com/samtools/samtools/releases/download/1.3/samtools-1.3.tar.bz2 && \
	tar -jxf samtools-1.3.tar.bz2 && \
	cd samtools-1.3 && \
	./configure --without-curses && \
	make && \
	make install	
#perl installation not neccessary since Development Tools, which yum installed earlier, includes Perl v5.10.1.
#INSTALL Perl. Needed for knife. On sherlock v5.10.1, I'll grab the latest, however.
#RUN mkdir /srv/src/Perl && \
#		cd /srv/src/Perl && \
#		wget http://www.cpan.org/src/5.0/perl-5.22.1.tar.gz && \
#		tar -zxf perl-5.22.1.tar.gz && \
#		cd perl-5.22.1 && \
#		sh Configure -de && \
#		make && 
#		make install &&

RUN git clone https://github.com/ericff/eknife.git /srv/software/knife
RUN git clone https://github.com/ericff/MACHETE.git /srv/software/machete

ENTRYPOINT []
LABEL version="1.0" description="Detects gene fusions"
