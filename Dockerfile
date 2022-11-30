FROM ubuntu:22.04 
RUN apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y \
      python3-pip \
      python3-venv \
      python3-astroquery \
      python3-astroplan \
      python3-future \
      python3-jinja2 \
      python3-matplotlib \
      python3-pandas \
      python3-patsy \
      python3-reportlab \
      python3-reproject \
      python3-seaborn \
      python3-threadpoolctl \
      python3-tornado \
      python3-tqdm \
    && rm -rf /var/cache/apt/* /var/lib/apt/lists/*
RUN useradd -m app
WORKDIR /home/app
COPY requirements.txt requirements.txt
USER app
RUN pip3 install --user -r requirements.txt
COPY . .
EXPOSE 5002
CMD ["/home/app/.local/bin/bokeh", "serve", ".", "--address=0.0.0.0", "--port=5002"]
