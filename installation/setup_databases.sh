#!/usr/bin/env bash
###############################################################################
# installation/setup_databases.sh
#
# Downloads/builds the reference databases that are too large to bundle
# directly in the GitHub repo, into databases/ at the repo root.
#
#   bash installation/setup_databases.sh
#
# NOT covered by this script — these are small, permissively-licensed, and
# already committed directly under databases/ in this repo, so nothing needs
# to be downloaded or built for them:
#   - databases/bacmet/            (BacMet2_EXP fasta + mapping file)
#   - databases/cog_funclass.tab
#   - databases/trimmomatic_adapters/
#
# Run individual sections by setting SKIP_<NAME>=1, e.g.:
#   SKIP_KAIJU=1 SKIP_EGGNOG=1 bash installation/setup_databases.sh
#
# PREREQUISITE: the `metores` conda environment must already exist and be
# active (see ../environment.yml) for wget/tar to be available as expected.
#
# SOURCES. Two official installers (Kaiju, eggNOG) are broken upstream;
# the per-section comments below say what is used instead:
#   Kaiju      : `kaiju-makedb -s progenomes` has a dead internal source URL
#                (github.com/bioinformatics-centre/kaiju/issues/273). This
#                script downloads the official PRE-BUILT index instead (AWS
#                S3, Open Data Sponsorship Program).
#   CARD       : https://card.mcmaster.ca/latest/data — CARD's stable
#                "always latest" download endpoint. NOTE: CARD's terms restrict *commercial*
#                reproduction of the data; research/non-commercial use is
#                unrestricted. See https://card.mcmaster.ca/about. This is
#                also why CARD is fetched here rather than bundled in git —
#                unlike BacMet, its license doesn't clearly permit that.
#   eggNOG     : the packaged `download_eggnog_data.py` hardcodes a dead
#                domain (eggnogdb.embl.de; eggnogdb/eggnog-mapper issues
#                #578, #589, #600). This script fetches the same files from
#                the current domain instead.
#   SortMeRNA  : the SortMeRNA v4.3.4 release tarball.
###############################################################################
set -uo pipefail  # NOT -e: one database failing shouldn't abort the rest

# Resolve the repo root regardless of where this script is invoked from
# (this file lives in installation/, so the repo root is one level up).
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
DB_DIR="${BASE_DIR}/databases"

mkdir -p "${DB_DIR}"

echo "=============================================="
echo " MeToReS — Database Setup"
echo " Repo root : ${BASE_DIR}"
echo " Databases : ${DB_DIR}"
echo "=============================================="

FAILED=""

# =============================================================================
# 1. Kaiju index (default: proGenomes, prokaryotes only)
# =============================================================================
# To remove plant/host reads you need an index that contains eukaryotes, e.g.
# nr_euk. Pick the URL from https://bioinformatics-centre.github.io/kaiju/downloads.html
# and run, for example:
#   KAIJU_INDEX=nr_euk KAIJU_URL=<url of the nr_euk .tgz> bash installation/setup_databases.sh
# then point kaiju.db / nodes / names in config.yaml at databases/kaiju/nr_euk/
# and set kaiju.taxon_filter (e.g. "Streptophyta"). Check the index's RAM
# requirement on that page and set kaiju.mem_mb accordingly.
if [ -z "${SKIP_KAIJU:-}" ]; then
    KAIJU_INDEX="${KAIJU_INDEX:-progenomes}"
    echo ""
    echo "[1/4] Kaiju (${KAIJU_INDEX})..."
    if [ "${KAIJU_INDEX}" = "progenomes" ]; then
        # kaiju-makedb's own proGenomes download step is broken
        # (github.com/bioinformatics-centre/kaiju/issues/273), so the
        # official pre-built index is downloaded instead. 2023-05-25 build:
        # ~66 GB compressed, ~102 GB RAM needed to run Kaiju with it.
        KAIJU_URL="${KAIJU_URL:-https://kaiju-idx.s3.eu-central-1.amazonaws.com/2023/kaiju_db_progenomes_2023-05-25.tgz}"
    elif [ -z "${KAIJU_URL:-}" ]; then
        echo "      KAIJU_INDEX=${KAIJU_INDEX} needs KAIJU_URL (see the downloads page above)."
        KAIJU_URL=""
        FAILED="${FAILED} kaiju"
    fi
    if [ -n "${KAIJU_URL}" ]; then
        mkdir -p "${DB_DIR}/kaiju/${KAIJU_INDEX}"
        (
            cd "${DB_DIR}/kaiju/${KAIJU_INDEX}" && \
            wget -q -O kaiju_index.tgz "${KAIJU_URL}" && \
            tar -xzf kaiju_index.tgz && \
            rm -f kaiju_index.tgz
        ) && echo "      Downloaded and extracted." \
          || { echo "      FAILED — check network access, or verify the URL is"; \
               echo "               still current at the Kaiju downloads page:"; \
               echo "               ${KAIJU_URL}"; \
               FAILED="${FAILED} kaiju"; }
        echo "      Expect kaiju_db_*.fmi, nodes.dmp, names.dmp in"
        echo "      ${DB_DIR}/kaiju/${KAIJU_INDEX}/"
    fi
else
    echo "[1/4] Kaiju: skipped (SKIP_KAIJU set)"
fi

# =============================================================================
# 2. CARD (for RGI / ARG annotation)
# =============================================================================
if [ -z "${SKIP_CARD:-}" ]; then
    echo ""
    echo "[2/4] CARD..."
    mkdir -p "${DB_DIR}/card"
    (
        cd "${DB_DIR}/card" && \
        wget -q -O card-data.tar.bz2 https://card.mcmaster.ca/latest/data && \
        tar -xjf card-data.tar.bz2 && \
        rm -f card-data.tar.bz2
    ) && echo "      card.json ready in ${DB_DIR}/card/" \
      || { echo "      FAILED — check network access to card.mcmaster.ca"; FAILED="${FAILED} card"; }
else
    echo "[2/4] CARD: skipped (SKIP_CARD set)"
fi

# =============================================================================
# 3. eggNOG-mapper database
# =============================================================================
if [ -z "${SKIP_EGGNOG:-}" ]; then
    echo ""
    echo "[3/4] eggNOG-mapper database (this is large, budget real time)..."
    mkdir -p "${DB_DIR}/eggnog-mapper"
    # NOTE: the packaged download_eggnog_data.py hardcodes the dead domain
    # eggnogdb.embl.de and fails with DNS/404 errors — this is a widely
    # reported, still-open upstream bug (e.g. github.com/eggnogdb/
    # eggnog-mapper issues #578, #589, #600), not specific to this server.
    # A fix exists as an unmerged PR. Rather than depend on that landing,
    # this fetches the same three files directly from the CORRECTED current
    # domain (eggnog5.embl.de), bypassing the broken script entirely.
    EGGNOG_BASE="http://eggnog5.embl.de/download/emapperdb-5.0.2"
    EGGNOG_OK=1
    (
        cd "${DB_DIR}/eggnog-mapper" && \
        wget -q -O eggnog.db.gz "${EGGNOG_BASE}/eggnog.db.gz" && \
        gunzip -f eggnog.db.gz
    ) || EGGNOG_OK=0
    (
        cd "${DB_DIR}/eggnog-mapper" && \
        wget -q -O eggnog.taxa.tar.gz "${EGGNOG_BASE}/eggnog.taxa.tar.gz" && \
        tar -zxf eggnog.taxa.tar.gz && \
        rm -f eggnog.taxa.tar.gz
    ) || EGGNOG_OK=0
    (
        cd "${DB_DIR}/eggnog-mapper" && \
        wget -q -O eggnog_proteins.dmnd.gz "${EGGNOG_BASE}/eggnog_proteins.dmnd.gz" && \
        gunzip -f eggnog_proteins.dmnd.gz
    ) || EGGNOG_OK=0

    if [ "${EGGNOG_OK}" -eq 1 ]; then
        echo "      Downloaded eggnog.db, eggnog.taxa.db, eggnog_proteins.dmnd."
    else
        echo "      FAILED — one or more files did not download. Check network"
        echo "               access to eggnog5.embl.de, or verify that domain"
        echo "               is still current (it has moved before)."
        FAILED="${FAILED} eggnog"
    fi
else
    echo "[3/4] eggNOG: skipped (SKIP_EGGNOG set)"
fi

# =============================================================================
# 4. SortMeRNA database
# =============================================================================
if [ -z "${SKIP_SORTMERNA:-}" ]; then
    echo ""
    echo "[4/4] SortMeRNA database..."
    mkdir -p "${DB_DIR}/sortmerna" "${DB_DIR}/sortmerna/idx"
    DB_FASTA="${DB_DIR}/sortmerna/smr_v4.3_default_db.fasta"
    if [ -s "${DB_FASTA}" ]; then
        echo "      Already present — skipping download."
    else
        (
            cd "${DB_DIR}/sortmerna" && \
            wget -q -O database.tar.gz \
                "https://github.com/sortmerna/sortmerna/releases/download/v4.3.4/database.tar.gz" && \
            tar -xzf database.tar.gz && \
            rm -f database.tar.gz
        ) && echo "      Downloaded and extracted." \
          || { echo "      FAILED"; FAILED="${FAILED} sortmerna"; }
    fi
    if [ -s "${DB_FASTA}" ]; then
        echo "      NOTE: this database is RNA alphabet (U, not T) — that is"
        echo "      EXPECTED for SortMeRNA and is handled internally. Do not"
        echo "      convert U to T (a DNA aligner given this file would fail"
        echo "      silently at 0% alignment; SortMeRNA itself handles it"
        echo "      correctly)."
    fi
else
    echo "[4/4] SortMeRNA: skipped (SKIP_SORTMERNA set)"
fi

# =============================================================================
# Summary
# =============================================================================
echo ""
echo "=============================================="
if [ -z "${FAILED}" ]; then
    echo " All databases set up successfully."
    echo " (BacMet, COG table, and Trimmomatic adapters were already present"
    echo " in databases/ — bundled with the repo, nothing to do for those.)"
else
    echo " Completed with issues in:${FAILED}"
    echo " Re-run individual sections after fixing the cause, e.g.:"
    echo "   SKIP_KAIJU=1 SKIP_CARD=1 bash installation/setup_databases.sh"
    echo "     # retries only eggNOG + SortMeRNA"
fi
echo ""
echo " Contents of ${DB_DIR}:"
ls -la "${DB_DIR}"
echo "=============================================="
