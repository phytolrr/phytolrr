FROM python:slim
COPY ./ /nsites/
ARG PIP_MIRROR
RUN pip install ${PIP_MIRROR} --no-cache-dir -r /nsites/requirements.txt
ENV PYTHONPATH /nsites
CMD python /nsites/web_service/lrr_search_web_service.py
