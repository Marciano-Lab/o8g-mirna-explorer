# o8G-miRNA Retargeting Explorer — container image
FROM python:3.13-slim

WORKDIR /app

# system deps kept minimal; wheels cover numpy/scipy/pyarrow
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# app code + precomputed data
COPY . .

# unpack gene-sets if only the tarball shipped
RUN if [ ! -d genesets ] && [ -f genesets.tar.gz ]; then tar -xzf genesets.tar.gz; fi

# GitHub rejects files over 100 MB. When an extract is shipped gzipped, expand it.
RUN if [ -f paper/data/hosted_refsets.sqlite.gz ] && [ ! -f paper/data/hosted_refsets.sqlite ]; then gzip -dc paper/data/hosted_refsets.sqlite.gz > paper/data/hosted_refsets.sqlite; fi
RUN if [ -f paper/data/Conserved_Site_Context_Scores.txt.gz ] && [ ! -f paper/data/Conserved_Site_Context_Scores.txt ]; then gzip -dc paper/data/Conserved_Site_Context_Scores.txt.gz > paper/data/Conserved_Site_Context_Scores.txt; fi
RUN if [ -f o8g_reverse.compact.db.gz ] && [ ! -f o8g_reverse.compact.db ]; then gzip -dc o8g_reverse.compact.db.gz > o8g_reverse.compact.db; fi

# managed platforms inject $PORT; default to 8501 locally
ENV PORT=8501
EXPOSE 8501

# 0.0.0.0 so the container is reachable; headless off the browser auto-open
CMD streamlit run app.py \
    --server.port=${PORT} \
    --server.address=0.0.0.0 \
    --server.headless=true \
    --server.enableCORS=false \
    --server.enableXsrfProtection=false
