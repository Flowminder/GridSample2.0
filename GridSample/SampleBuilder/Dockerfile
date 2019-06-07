##
# Dockerfile for grid sample builder
#
# This is based on the geodata/gdal container
# and adds python,postgres/gis and the python modules required to run the 
# Grid Sample python code 

FROM geodata/gdal:2.1.2

MAINTAINER Ian Waldock <icw@geodata.soton.ac.uk>

# Install the dependencies.
# r-base r-base-core : to run R programs - Sample builder runs r scripts from within it
# although we  build r from source later install them anyway as it makes building from source work ok as it seems to install dependencies such as a fortran compiler!
USER root
RUN apt-get update -y && \
    apt-get install -y sshfs \
      postgresql-client \
      postgis \
      python-psycopg2 \
      p7zip-full \
	cron \
      wget \
      python-pip \
	python-numpy python-scipy python-matplotlib ipython ipython-notebook python-pandas python-sympy python-nose \
	sendmail r-base r-base-core wkhtmltopdf


#RUN python -m pip install geographiclib
#RUN python -m pip install shapely
#RUN python -m pip install fiona
#RUN python -m pip install six==1.10.0


# install r from source to get latest version - as some required packages need a latest version
RUN wget http://www.stats.bris.ac.uk/R/src/base/R-3/R-3.5.2.tar.gz \
&& tar -xvzf R-3.5.2.tar.gz \
&& cd R-3.5.2 \
&& ./configure --enable-R-shlib \
&& make \
&& make install 

# install r packages that we need for our app - the 
RUN R -e "install.packages('sp', method='wget',repos='http://www.stats.bris.ac.uk/R/')"
RUN R -e "install.packages('rgdal', method='wget',repos='http://www.stats.bris.ac.uk/R/')"
RUN R -e "install.packages('raster', method='wget',repos='http://www.stats.bris.ac.uk/R/')"
RUN R -e "install.packages('igraph', method='wget',repos='http://www.stats.bris.ac.uk/R/')"

# modules to run r scripts from python code
RUN python -m pip install singledispatch
RUN python -m pip install rpy2==2.5.2
#RUN python -m pip install rpy2

# to create pdf from html
RUN python -m pip install pdfkit

# virtual x-term for html to pdf lib
RUN apt-get install -y xvfb \
&& printf '#!/bin/bash\nxvfb-run -a --server-args="-screen 0, 1024x768x24" /usr/bin/wkhtmltopdf -q $*' > /usr/bin/wkhtmltopdf.sh \
&& chmod a+x /usr/bin/wkhtmltopdf.sh \
&& ln -s /usr/bin/wkhtmltopdf.sh /usr/local/bin/wkhtmltopdf

# we need to chain the commands as there seems to be no permanent dir structure when building
RUN wget https://www.libarchive.org/downloads/libarchive-3.3.1.tar.gz \
&& tar xzf libarchive-3.3.1.tar.gz \
&& cd libarchive-3.3.1 \
&& ./configure \
&& make \
&& make install \
&& ldconfig 

# to install libarchive to manipulate .7z file
RUN python -m pip install libarchive

#RUN cd ..
#RUN rm -rf libarchive-3.3.1
#RUN rm libarchive-3.3.1.tar.gz


# set up crontab
# Add crontab file in the cron directory
#ADD src/cron.txt /cron.txt
#RUN crontab /cron.txt
#RUN rm /cron.txt

# set up log file
#RUN touch /var/log/cron.log

#COPY ./src /src
#WORKDIR /src

# set lib path so python can call R scripts
ENV LD_LIBRARY_PATH /usr/local/lib/R/library/stats/libs/:/usr/local/lib/R/lib/

# just to keep the container alive
CMD tail -f /dev/null
#CMD ["/src/startUp.sh"]