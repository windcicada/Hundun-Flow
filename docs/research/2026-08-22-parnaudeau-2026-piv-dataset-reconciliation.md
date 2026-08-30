# Parnaudeau 2026 PIV supplementary-dataset reconciliation

**Status:** NOT_EQUIVALENT for replacing the frozen Parnaudeau Fig. 11--15
controlled digitization. This is a read-only research receipt only. It does not
modify any literature JSON, receipt completion flag, plan, ledger, solver,
COAST lane, or Stage 5 artifact, and it does not authorize a tolerance,
full20, long statistics, or a final gate.

**Scope and freeze:** The frozen reconciliation instruction is
docs/verification/v0.4-literature-data-receipt.md at commit c11c73d. The only
dataset inspected is Recherche Data Gouv persistent ID
doi:10.57745/DHJXM6, released V1.0. The official API inventory was enumerated
before any file data was fetched. The approximately 12.3 GB instantaneous
archive was not downloaded. No HUNDUN benchmark, output, or evidence was read.

## Result matrix

The statuses below are deliberately conservative. A positive statement in a
source is recorded as fact; a missing or ambiguous field is not filled from a
secondary CFD source.

| frozen row | status | direct official locator and short source text | evidence and consequence |
|---|---|---|---|
| experiment identity | UNRESOLVED | Dataverse API /data/latestVersion/metadataBlocks/citation/fields, dsDescriptionValue: “experimental data ... circular cylinder at Reynolds number 3900”; “2D2C PIV ... in an INRAE wind tunnel”; and “completes and improves ... first and second order statistics ... provided in Parnaudeau et al. (2008)” | The released metadata explicitly links this dataset to the Parnaudeau 2008 PIV database. It does not supply cylinder diameter, reference velocity, blockage, tunnel dimensions, or the PIV plane. Thus the relation is a documented claim, but the complete experiment identity needed for substitution remains unresolved. |
| coordinate identity | NOT_EQUIVALENT | First released file, byte range 0--8191, first line: #DaVis 10.x 2C vector field 4 545 740 "x [mm]";"y [mm]";"Vx [m/s]";"Vy [m/s]" | The header gives x/y labels and millimetres, but no physical origin, cylinder-centre convention, diameter, axis orientation relative to the paper, or explicit mapping to x/D=1.06,1.54,2.02. A station mapping without extrapolation cannot be established. |
| normalization identity | NOT_EQUIVALENT | The same authoritative header labels Vx and Vy as [m/s]; dataset metadata gives Re=3900 but no Uc or normalization rule | Physical velocity units are present in a raw field header, but the source value and definition required for conversion to Uc and Uc^2 are absent. No normalized comparison may be inferred. |
| observable coverage | NOT_EQUIVALENT | Dataverse /data/latestVersion/files contains 998 entries, all contentType=text/plain, all filename-matching Serie_*.txt; no published mean/stress product is enumerated | No released machine-readable product was found that supplies all five required quantities (mean u, mean v, u'u', v'v', u'v') at all three registered stations over a common transverse interval. The descriptive metadata claims that mean and Reynolds-stress fields exist, but the V1 file inventory exposes only the Serie files; the claim cannot substitute for products. |
| extraction error | UNRESOLVED | No mean/Reynolds-stress file, README, station table, or extraction protocol is present in the released V1 inventory; header has a finite grid declaration but no station mapping | Exact-grid selection, interpolation rule, station offset, and extraction/interpolation error cannot be computed or pre-registered. No final tolerance is designed here. |
| provenance | UNRESOLVED | Dataverse V1 metadata, file endpoint, DataCite DOI record, and the local Range artifact are listed below | Persistent ID, version, release time, official file IDs/names/sizes/MD5 values, and license are bound. Full-file local SHA256 values are intentionally unavailable because the raw archive was not downloaded; only the 8192-byte header Range has a local SHA256. Provenance is therefore not sealed for a profile replacement. |

**Overall decision:** NOT_EQUIVALENT. The explicit Parnaudeau link is useful
provenance, but the missing coordinate, normalization, published five-observable
coverage, and local full-file evidence prevent this dataset from replacing the
frozen controlled digitization.

## Official dataset identity and version

- Persistent ID: doi:10.57745/DHJXM6.
- Landing page: https://entrepot.recherche.data.gouv.fr/dataset.xhtml?persistentId=doi%3A10.57745%2FDHJXM6
- Persistent URL: https://doi.org/10.57745/DHJXM6
- Official Dataverse publisher: Recherche Data Gouv.
- V1.0, state RELEASED; dataset version ID 283288; dataset ID 720391;
  release/publication time 2026-02-19T18:24:44Z.
- Dataverse license: etalab 2.0,
  https://spdx.org/licenses/etalab-2.0.html.
- Dataset title: “Non-time-resolved PIV dataset of flow over a circular
  cylinder at Reynolds number 3900”.
- Official contact in citation metadata: Dominique Heitz, INRAE,
  dominique.heitz@inrae.fr.
- DataCite confirms DOI 10.57745/dhjxm6, publisher Recherche Data Gouv,
  publication year 2026, and 998 HasPart child DOIs. DataCite does not provide
  a version or license field for this record; the license above is from the
  official Dataverse V1 record.

### Published description versus released file inventory

The official description says the dataset contains 10,000 instantaneous
velocity fields, mean velocity and Reynolds-stress fields, and skewness and
flatness fields, and says it completes/improves the Parnaudeau 2008
first/second-order PIV database. The same V1 API complete files array has 998
entries, all text/plain, with names only of the form Serie_NNNNN.txt; no README
or mean/stress/skewness/flatness filename or file-level description occurs. The
sequence is 10001--10999 with 10902 absent. The inventory therefore records
what is actually downloadable in V1, without trying to reconcile the
descriptive count by guessing.

| inventory fact | official value |
|---|---|
| API file count | 998 |
| filename pattern | Serie_10001.txt through Serie_10999.txt, except Serie_10902.txt |
| extension/content type | 998 .txt / text/plain |
| total advertised bytes | 13,176,725,066 |
| advertised size range | 13,156,978--13,252,930 bytes |
| child DOI count in DataCite | 998 HasPart DOIs |
| small published profile/readme product | none found |
| archive action | not downloaded; only one 8192-byte HTTP Range header fetched |

The full official inventory (ID, filename, advertised byte size, official MD5,
child persistent ID, and datafile endpoint) is reproduced in the appendix. Its
TSV source is also retained for this session at
/tmp/dhjxm6-file-inventory.tsv, SHA256
07b12407247ffc414383adcc0092b2d4d020015eab07f78c15115e0641380044.

Representative official file records (the appendix is authoritative for all
998 entries):

| file ID | filename | size bytes | official checksum | child persistent ID |
|---:|---|---:|---|---|
| 720393 | Serie_10001.txt | 13,198,574 | MD5 ad8e65c5f4e028763db96a957629d644 | doi:10.57745/LPNPL5 |
| 720401 | Serie_10002.txt | 13,216,263 | MD5 0a1fe9a93bf1d785e3e84ce5c53ecf03 | doi:10.57745/ATEWCQ |
| 720399 | Serie_10003.txt | 13,199,721 | MD5 2ec95948db3f00ab4c8fe760b336d0c0 | doi:10.57745/QCXFQH |
| 720929 | Serie_10900.txt | 13,174,138 | MD5 4c617190b4305d125fd15db1771f3ce5 | doi:10.57745/OCFPIW |
| 721154 | Serie_10901.txt | 13,187,532 | MD5 5788bb7563b3c8b3fcf6948ad5a11ae0 | doi:10.57745/XERZCU |
| — | Serie_10902.txt | — | — | missing from V1 inventory |
| 721370 | Serie_10903.txt | 13,214,274 | MD5 8740531a64b53769e94ffa6c5b9be124 | doi:10.57745/E0I16U |
| 721029 | Serie_10998.txt | 13,208,427 | MD5 418446757a16b4f2044145694b7eea0e | doi:10.57745/P3QLXY |
| 721086 | Serie_10999.txt | 13,216,581 | MD5 34cfb2ff14efdf2cd049a2ec4fc766c5 | doi:10.57745/0KEKUM |

## Header, coordinates, and units

Only an HTTP Range of the first file was downloaded:

- Official endpoint:
  https://entrepot.recherche.data.gouv.fr/api/access/datafile/720393.
- Request: Range: bytes=0-8191; local file
  /tmp/dhjxm6-serie10001-header.bin.
- Download mtime: 2026-08-22 06:43:48.942923994 +0800.
- Local artifact size: 8192 bytes.
- Local artifact SHA256:
  9b3b5974781080bdaedb160a5c5c07b9a52a06853ec90c1cd99ff65392a2350b.
- Official whole-file record: ID 720393, Serie_10001.txt, 13,198,574
  bytes, MD5 ad8e65c5f4e028763db96a957629d644.
- Exact first header line (short factual excerpt):
  #DaVis 10.x 2C vector field 4 545 740 "x [mm]";"y [mm]";"Vx [m/s]";"Vy [m/s]".
- The first displayed coordinate row is -60.1328;60.1803;0;-0.

The header establishes a 4-component DaVis vector-field declaration with a
545-by-740 grid and x/y in millimetres, Vx/Vy in metres per second. It does
not establish a physical origin, cylinder-centre location, D, streamwise
direction convention, PIV plane, or the conversion from a header coordinate to
the paper's x/D stations.

## Experiment parameters

| parameter required by frozen reconciliation | official source value | status |
|---|---|---|
| cylinder | circular cylinder, from title/description | identity link present; diameter absent |
| Reynolds number | 3900 | supplied |
| reference velocity | not present in V1 API/landing metadata or sampled header | UNRESOLVED |
| cylinder diameter D | not present | UNRESOLVED |
| blockage ratio | not present | UNRESOLVED |
| wind-tunnel dimensions/geometry | only “INRAE wind tunnel” is stated | UNRESOLVED |
| PIV plane | not present | UNRESOLVED |
| method | planar 2D2C PIV, non-time-resolved | supplied |
| normalization | raw Vx/Vy units m/s; no Uc definition/value | NOT_EQUIVALENT |

The primary paper is bound only as primary bibliographic/abstract provenance,
not as a source for unreported dataset parameters:

- HAL record: https://hal.science/hal-00383669v1 and official HAL API
  https://api.archives-ouvertes.fr/search/?q=halId_s%3Ahal-00383669&wt=json.
- DOI: https://doi.org/10.1063/1.2957018, Physics of Fluids 20 (2008),
  article 085101, “Experimental and numerical studies of the flow over a
  circular cylinder at Reynolds number 3900”.
- HAL doiId_s and en_abstract_s identify the paper and state that it
  investigates Re=3900 with LES, hot-wire anemometry, and PIV. The abstract
  does not provide the missing D, Uref, blockage, tunnel geometry, PIV-plane
  origin, or the three station arrays. No secondary CFD values were used.

## Five observables at three registered stations

The required machine-readable observables are:
u, v, u'u', v'v', and u'v', each at x/D=1.06, 1.54, 2.02, over a common
transverse interval.

| observable | x/D=1.06 | x/D=1.54 | x/D=2.02 | source/effect |
|---|---|---|---|---|
| mean u | NOT_EQUIVALENT | NOT_EQUIVALENT | NOT_EQUIVALENT | no released mean-profile file |
| mean v | NOT_EQUIVALENT | NOT_EQUIVALENT | NOT_EQUIVALENT | no released mean-profile file |
| u'u' | NOT_EQUIVALENT | NOT_EQUIVALENT | NOT_EQUIVALENT | no released Reynolds-stress file |
| v'v' | NOT_EQUIVALENT | NOT_EQUIVALENT | NOT_EQUIVALENT | no released Reynolds-stress file |
| u'v' | NOT_EQUIVALENT | NOT_EQUIVALENT | NOT_EQUIVALENT | no released Reynolds-stress file |

Common transverse coverage, numerical arrays, station pairing, and any
normalization to Uc/Uc^2 are therefore absent. The descriptive metadata's
assertion that mean/stress fields exist is not machine-readable coverage and is
not promoted.

## Station bias and extraction error

- Station targets: x/D=1.06, 1.54, 2.02.
- Dataset-grid station coordinates: not computable; D and the physical origin
  are absent.
- Maximum station offset: unresolved, not zero by assumption.
- Exact-grid selection/interpolation rule: not selected because the station
  mapping is not defined.
- Resulting extraction/interpolation error: unresolved.
- Experimental uncertainty: not inferred from the 2026 metadata, and not
  conflated with the Parnaudeau paper's abstract statement about statistical
  estimation.
- No comparison tolerance is designed in this receipt. A future route would
  first need an official profile product or a separately frozen raw-data
  reduction protocol, then a pre-registered exact-grid/interpolation rule and
  an independently calculated extraction error.

## Provenance and downloaded artifacts

### Official source snapshots

| artifact | official URL | local path | local SHA256 | download mtime (+08:00) |
|---|---|---|---|---|
| Dataverse dataset V1 API | https://entrepot.recherche.data.gouv.fr/api/datasets/:persistentId/?persistentId=doi%3A10.57745%2FDHJXM6 | /tmp/dhjxm6-dataset-api.json | 02f0a612e515788073e3380be99f21df73e18a01634fbd81afaecaa9c2f55aed | 2026-08-22 06:40:45.379062345 |
| Dataverse versions API | https://entrepot.recherche.data.gouv.fr/api/datasets/:persistentId/versions?persistentId=doi%3A10.57745%2FDHJXM6 | /tmp/dhjxm6-versions.json | f74a7b1eb6b4bde084cc30dfbbc4f806ccd8c4367668a7068aca72ba025b48a0 | 2026-08-22 06:45:27.894887489 |
| Dataverse file metadata (720393) | https://entrepot.recherche.data.gouv.fr/api/files/720393 | /tmp/dhjxm6-file-720393.json | 377f3a942ad4d0fe13260e4260b962e4ea33498fc644ccf771b231014e0382a3 | 2026-08-22 06:45:35.215164163 |
| official landing page | https://entrepot.recherche.data.gouv.fr/dataset.xhtml?persistentId=doi%3A10.57745%2FDHJXM6 | /tmp/dhjxm6-landing.html | 28c307b8da5d300f534309ab890670844c267385c527597bd3ea5f9c83aa09ea | 2026-08-22 06:41:55.202052160 |
| DataCite DOI API | https://api.datacite.org/dois/10.57745/DHJXM6 | /tmp/dhjxm6-datacite.json | ec38a641e21a1e44adb6f17cdfca33f0bb7cdb2c8da08e5822ce28511f4c5184 | 2026-08-22 06:44:16.704113331 |
| primary-paper HAL API | https://api.archives-ouvertes.fr/search/?q=halId_s:hal-00383669&wt=json&fl=* | /tmp/parnaudeau-hal.json | 16ddc8efffe30eba4bc98873d392bb64a9e072b407f2fd953cac2515077a58a9 | 2026-08-22 06:43:07.641154732 |
| primary-paper HAL file-field query | https://api.archives-ouvertes.fr/search/?q=halId_s:hal-00383669&wt=json&fl=uri_s,file*,docType_s,*Url* | /tmp/parnaudeau-hal-files.json | 85edc58957b0bdb7eb578e3b3e0271fd58e6484e29ccc54eba4bd540ba516b72 | 2026-08-22 06:43:19.193649588 |
| DOI bibliographic API | https://api.crossref.org/works?query.title=Parnaudeau%20flow%20around%20a%20circular%20cylinder%20at%20Reynolds%20number%203900&rows=5 | /tmp/parnaudeau-crossref.json | 72b1a37ce8bcd8b044d331d442d37786c84ee11f8d9d852034e0842515aa09b4 | 2026-08-22 06:42:17.391002447 |

The AIP PDF endpoint returned access-denied in this environment; the official
HAL metadata was sufficient to bind the primary paper's DOI/title/abstract.

### Temporary-file inventory

- /tmp/dhjxm6-dataset-api.json — complete V1 API response, 499,910 bytes,
  SHA256 above.
- /tmp/dhjxm6-versions.json — V1 versions response, 499,635 bytes, SHA256
  above.
- /tmp/dhjxm6-file-720393.json — authoritative file metadata for sampled
  header, SHA256 above.
- /tmp/dhjxm6-landing.html — official landing page, 516,652 bytes, SHA256
  above.
- /tmp/dhjxm6-datacite.json — DataCite DOI record, SHA256 above.
- /tmp/parnaudeau-hal.json and /tmp/parnaudeau-hal-files.json — official HAL
  metadata queries, SHA256 above.
- /tmp/dhjxm6-serie10001-header.bin — 8192-byte Range artifact, SHA256 above.
- /tmp/dhjxm6-file-inventory.tsv — all 998 official file records, 151,750
  bytes, SHA256
  07b12407247ffc414383adcc0092b2d4d020015eab07f78c15115e0641380044.
- No instantaneous archive or full Serie file was downloaded.

## Reproduction commands and tool hashes

All commands were read-only HTTP GET/range requests or local parsing/hashing.
The working directory was /home/wyf/code_dev/hundun-flow; no benchmark,
output, evidence, or solver path was accessed.

~~~sh
api='https://entrepot.recherche.data.gouv.fr/api/datasets/:persistentId/?persistentId=doi%3A10.57745%2FDHJXM6'
curl -fsSL --retry 2 --connect-timeout 20 "$api" -o /tmp/dhjxm6-dataset-api.json
curl -fsSL --retry 2 --connect-timeout 20 'https://entrepot.recherche.data.gouv.fr/dataset.xhtml?persistentId=doi%3A10.57745%2FDHJXM6' -o /tmp/dhjxm6-landing.html
curl -fsSL --retry 2 --connect-timeout 20 'https://api.datacite.org/dois/10.57745/DHJXM6' -o /tmp/dhjxm6-datacite.json
curl -fsSL --retry 2 --connect-timeout 20 'https://api.archives-ouvertes.fr/search/?q=halId_s:hal-00383669&wt=json&fl=*' -o /tmp/parnaudeau-hal.json
curl -fsSL --retry 2 --connect-timeout 20 -H 'Range: bytes=0-8191' \
  'https://entrepot.recherche.data.gouv.fr/api/access/datafile/720393' \
  -o /tmp/dhjxm6-serie10001-header.bin
python3 - <<'PY'
import json
j = json.load(open('/tmp/dhjxm6-dataset-api.json'))
fs = sorted((x['dataFile'] for x in j['data']['latestVersion']['files']),
            key=lambda f: f['filename'])
with open('/tmp/dhjxm6-file-inventory.tsv', 'w', encoding='utf-8') as out:
    out.write('id\\tfilename\\tsize_bytes\\tmd5\\tpersistent_id\\tdatafile_url\\n')
    for f in fs:
        out.write(f"{f['id']}\\t{f['filename']}\\t{f['filesize']}\\t"
                  f"{f['md5']}\\t{f['persistentId']}\\t"
                  f"https://entrepot.recherche.data.gouv.fr/api/access/datafile/{f['id']}\\n")
PY
sha256sum /tmp/dhjxm6-dataset-api.json /tmp/dhjxm6-landing.html \
  /tmp/dhjxm6-datacite.json /tmp/parnaudeau-hal.json \
  /tmp/parnaudeau-hal-files.json /tmp/dhjxm6-serie10001-header.bin \
  /tmp/dhjxm6-file-inventory.tsv
~~~

Relevant tool versions and executable SHA256:

| tool | version | executable SHA256 |
|---|---|---|
| curl | 7.68.0 | c1e06a1d2e28f984f9d4d338560423f9c78c4cddccafd3e0c28edfa55fcb708b |
| Python | 3.6.9 | 5e0c61520f859fea220cbe491c687b7a41f380050a5bcfd78b5917bcec921814 |
| sha256sum | GNU coreutils | 08ad66a3d429596f28c674074b6bbedaddb3025a16275ddcd72a7d6058bd5a80 |
| sed | GNU | 8d109b799969593c9453226e8a5ee8ba1762bc049a366771fcaa0e2365699c08 |
| od | GNU | 126042c17073409307854500893c35372f9310b113c0b6497148368934faae22 |
| file | 5.32 | 96ccb3180c1b62a9bb6910d7505190e15716ed109279a6d036ba020377018958 |
| ripgrep | 15.2.0 | e62198eb19b136b88c330af83647b5a962cb99b6b1f066758568f12de1974849 |

## Appendix: complete V1 file inventory

The following 998 rows are copied from the successful official Dataverse API
response and are not re-derived from any HUNDUN or secondary CFD artifact.

~~~tsv
id	filename	size_bytes	md5	persistent_id	datafile_url
720393	Serie_10001.txt	13198574	ad8e65c5f4e028763db96a957629d644	doi:10.57745/LPNPL5	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720393
720401	Serie_10002.txt	13216263	0a1fe9a93bf1d785e3e84ce5c53ecf03	doi:10.57745/ATEWCQ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720401
720399	Serie_10003.txt	13199721	2ec95948db3f00ab4c8fe760b336d0c0	doi:10.57745/QCXFQH	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720399
720392	Serie_10004.txt	13224196	3341714cea3eeb594a377e8b9ccb67ad	doi:10.57745/FY2FTH	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720392
720397	Serie_10005.txt	13234737	0174db0e435f16e67eee8e381e1e7faa	doi:10.57745/UOFATV	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720397
720400	Serie_10006.txt	13193442	c3a06c847fdf39b5094ebf860e350e2a	doi:10.57745/PWPLJQ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720400
720394	Serie_10007.txt	13212039	a8c30ed539d1caefd4aba8ad3c19c686	doi:10.57745/BAJTAF	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720394
720398	Serie_10008.txt	13214886	909da671ce13e6aaf0a006b3798db669	doi:10.57745/4ZOTT8	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720398
720395	Serie_10009.txt	13203598	d4d91f554f60699fabe18c9ffde395f8	doi:10.57745/SHV2TB	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720395
720396	Serie_10010.txt	13190466	b6b3395c2e2383f741ec158297e8e970	doi:10.57745/MY87GX	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720396
721322	Serie_10011.txt	13194046	da8e0409d398e70f53ae1abf6cba3bb9	doi:10.57745/VEOEQO	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721322
721244	Serie_10012.txt	13252930	d9d55961b5da8de82c7ea0a7d883aa57	doi:10.57745/Y7WXI5	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721244
721305	Serie_10013.txt	13196940	4b3a3319d6552ca4d98a33d5fd35b479	doi:10.57745/AOJH9C	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721305
720869	Serie_10014.txt	13179852	4dcec06a14953554d323ad590c665821	doi:10.57745/EOMDTX	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720869
720575	Serie_10015.txt	13189973	5011acecaa06aa7bd64995ee1cd05161	doi:10.57745/3RUI0D	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720575
720738	Serie_10016.txt	13212946	55f82113bf238a052fa120bc8155453d	doi:10.57745/8AYLBK	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720738
721299	Serie_10017.txt	13199884	e6301ea1c04e4d6331c0d001f482a530	doi:10.57745/FTR9RM	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721299
720630	Serie_10018.txt	13183960	8929e081b2e88377c2a27846da68842a	doi:10.57745/IEZ6IT	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720630
720693	Serie_10019.txt	13210959	731cc20cdb21ddbc65f01eb62fe1b954	doi:10.57745/ZYGEHO	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720693
720864	Serie_10020.txt	13197007	32de12374b1fb6f14ddf05067ef8bd9d	doi:10.57745/9DDI3Q	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720864
721111	Serie_10021.txt	13197258	e499aa66dbb388be64bc2dfde0f189df	doi:10.57745/ZQU9ZU	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721111
720851	Serie_10022.txt	13185278	e3d8675adf1c57553ba79718dc30eafc	doi:10.57745/TIGOWI	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720851
720975	Serie_10023.txt	13219432	85fa1c42a74866a4c1d9793990186ba3	doi:10.57745/6LDFFI	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720975
721152	Serie_10024.txt	13218257	eb78d441977d8458610272513604fd44	doi:10.57745/CZDOTM	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721152
721187	Serie_10025.txt	13185878	cd9d158f289f5bef8ae8de8d3460b7ec	doi:10.57745/7XS34U	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721187
721043	Serie_10026.txt	13217080	5088266b403539b8c2c0d236c821d38b	doi:10.57745/JLKMEO	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721043
720753	Serie_10027.txt	13217898	bcc51a987488c619930de89c64318e2c	doi:10.57745/PPQDL7	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720753
720934	Serie_10028.txt	13178072	3618eadc70cc48c56aa181416b25dd93	doi:10.57745/31932O	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720934
720802	Serie_10029.txt	13204945	47b35c1991c043de3021730bce7b15e0	doi:10.57745/RPAGM4	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720802
721023	Serie_10030.txt	13214709	46fb52cdad0b2a8e489500c953c49101	doi:10.57745/EDOHLH	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721023
721383	Serie_10031.txt	13168580	055662b9af6db698449f4e3bb09c2091	doi:10.57745/QLJLLS	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721383
720467	Serie_10032.txt	13227253	b88c3e5c8b1c2e190a5392053e52568d	doi:10.57745/RSCOYD	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720467
721004	Serie_10033.txt	13206296	4e126a11a9689f61fe95ef4c88624116	doi:10.57745/SS2REU	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721004
721238	Serie_10034.txt	13191635	46985db0a957188035d3786e8f7ac849	doi:10.57745/ONWE1M	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721238
720977	Serie_10035.txt	13201296	da8291f55eecae14049b2fbf8420f575	doi:10.57745/VIHBY8	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720977
720448	Serie_10036.txt	13237702	2b3160dae0cec5e72c6c42a5e323a21e	doi:10.57745/CQHQOX	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720448
720912	Serie_10037.txt	13229411	4d5049c18e8b76e2e8d50b629dcca176	doi:10.57745/S1YWDO	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720912
720905	Serie_10038.txt	13216142	befe6abc6149984d61913ce8bccb2e9e	doi:10.57745/FANQBC	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720905
721113	Serie_10039.txt	13206719	fcadff2c674159ce73c7fb8cbbdb4944	doi:10.57745/VMCWYV	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721113
720657	Serie_10040.txt	13177941	ae992f7f71c13591a7848ec14ff2fb2d	doi:10.57745/7JYLBP	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720657
720637	Serie_10041.txt	13195385	df1094f5fcb9dc275c47d7bfccb04863	doi:10.57745/WSDQSC	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720637
721020	Serie_10042.txt	13204273	ad34f134468bd81f7fcf7e04ae3d350e	doi:10.57745/WTRQTU	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721020
721330	Serie_10043.txt	13184705	f5a3f9cb27ef54863a092c8e6667345e	doi:10.57745/JRB2MU	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721330
720750	Serie_10044.txt	13188660	df6ad5ddf72d3c3a23fef34567824439	doi:10.57745/GHSSOO	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720750
721070	Serie_10045.txt	13206100	4006c874650aceeb18f4ecabf5f80704	doi:10.57745/YVTRPH	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721070
721328	Serie_10046.txt	13198111	d79b15edc0c44d87bf969c47a14d5d48	doi:10.57745/JULSKX	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721328
720760	Serie_10047.txt	13207758	6584b405b8ec100fe593c0440cb48a47	doi:10.57745/83BFT4	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720760
720572	Serie_10048.txt	13208340	b4dd8f56c11b48d1428576eca5c11fae	doi:10.57745/DDSMUN	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720572
720493	Serie_10049.txt	13187980	599c1455a375a0c8f38d2f6321cc673b	doi:10.57745/TFUKEX	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720493
721231	Serie_10050.txt	13187809	9c264cc8eaf51bc20b5ab65c2e20be47	doi:10.57745/IFKBM0	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721231
720464	Serie_10051.txt	13200795	8c49bebd122866639c416ecbd59bf2c2	doi:10.57745/O6VMGS	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720464
721037	Serie_10052.txt	13191251	95a541232eab4c9531d340e2c66b2a71	doi:10.57745/OZLGHR	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721037
720677	Serie_10053.txt	13200054	17ca2006ebae82f2f0286dbe4a53f9a9	doi:10.57745/4DDFYT	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720677
720528	Serie_10054.txt	13207643	655720de20a7706f8c8219bc2524d6ed	doi:10.57745/8XGXL5	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720528
721224	Serie_10055.txt	13213099	ca22eefb043fb38e2d901da58294d380	doi:10.57745/GD9CL1	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721224
720781	Serie_10056.txt	13183473	05a17a8b4f741a23432cc52f8a786d28	doi:10.57745/GMVTAT	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720781
721017	Serie_10057.txt	13186262	3ac45328007b66627a259e454a599164	doi:10.57745/ALQUNU	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721017
720612	Serie_10058.txt	13171316	5baf19328cf5dda72cc3c9bced23e586	doi:10.57745/NOAXGG	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720612
721250	Serie_10059.txt	13192660	78a2e0732456d62b69b9230d1ef5d398	doi:10.57745/PTLRCH	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721250
720950	Serie_10060.txt	13231791	d12d2014b1e2cbd0923afbdf49a968cf	doi:10.57745/BZP5LE	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720950
720730	Serie_10061.txt	13193943	a4b630a984edc7066e17737d1d96ef87	doi:10.57745/T1YV9L	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720730
720847	Serie_10062.txt	13195677	a241281b017867216b39839fd23b7232	doi:10.57745/V0MNVO	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720847
721040	Serie_10063.txt	13225753	adb52200a6d9d3372c36700ea931fc9d	doi:10.57745/DVOB0E	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721040
720849	Serie_10064.txt	13215101	c3a1c608b4bd8d47f87e319456cdd64e	doi:10.57745/ASKPM6	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720849
720999	Serie_10065.txt	13196364	094504aa29e26e610b0dbed3bec1d35d	doi:10.57745/NCQJOV	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720999
720902	Serie_10066.txt	13183096	299142385375ab2aac3f167e2d5f6c63	doi:10.57745/FBD0XN	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720902
720672	Serie_10067.txt	13174958	a6320e4e63ce2ecec12a896990c8e8ed	doi:10.57745/MSOAPJ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720672
720879	Serie_10068.txt	13199106	b6b6350081d18d57dfaf8a7c4a8f0458	doi:10.57745/GW933K	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720879
720811	Serie_10069.txt	13204841	872525b38887258fb546f07056fa5ef4	doi:10.57745/DWXDQX	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720811
720646	Serie_10070.txt	13218969	2872677e465e9d1dd43348383e29ebe2	doi:10.57745/I0ZPCB	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720646
720487	Serie_10071.txt	13201882	4e140d97513798c31db6c4115b86ab78	doi:10.57745/WVJIO7	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720487
721206	Serie_10072.txt	13219667	11f6c295afda160c0c1f52ffd7b6f41c	doi:10.57745/RU1XNZ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721206
721373	Serie_10073.txt	13178031	8918daadf108dd448a5eba95e338298d	doi:10.57745/T6FYKA	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721373
720697	Serie_10074.txt	13181324	f3c485c60b925b5bd27a3ab67bdf81bf	doi:10.57745/BAPZ39	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720697
720427	Serie_10075.txt	13195473	aa3978554899929912085a3874becdfd	doi:10.57745/5059CT	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720427
721237	Serie_10076.txt	13195183	0bebc951c082d91344c99ace6834683d	doi:10.57745/V6E8CP	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721237
721108	Serie_10077.txt	13197990	a9cfbb64957d24b41744744c75c11810	doi:10.57745/RGXBCX	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721108
720479	Serie_10078.txt	13195457	f3278d0e3fa91564817c2c04ad07f683	doi:10.57745/GWLTPJ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720479
721232	Serie_10079.txt	13221046	c139b372f02f80b0b5846bc705aa85c7	doi:10.57745/QJAJYD	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721232
720874	Serie_10080.txt	13211693	bb14fb976d75b6f1ff7f57f310ee3313	doi:10.57745/OUGHB9	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720874
720940	Serie_10081.txt	13218948	0de6abe5e71229754d7d8609fbfb7ddf	doi:10.57745/TYDTJM	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720940
720477	Serie_10082.txt	13199335	ca258f28908b7974e717c24d7aa5dfbe	doi:10.57745/UBQFVC	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720477
721302	Serie_10083.txt	13225688	b35e5e4bb314c511fb9a9b0e34dfe4be	doi:10.57745/CXBWSF	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721302
720446	Serie_10084.txt	13220458	29a1314936f47b6ae2cbb291c39aa059	doi:10.57745/KLI0SM	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720446
720696	Serie_10085.txt	13218389	c0dac453ee78879fcb21ad5b355a9e9b	doi:10.57745/65EN6E	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720696
720638	Serie_10086.txt	13221258	7fc0de61c9ddb3ade0e6d1f2a8439ea7	doi:10.57745/NDWXZL	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720638
720728	Serie_10087.txt	13222519	b1f2a8a50b7c7a525775a7861a4c2dd1	doi:10.57745/B7OFOT	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720728
720827	Serie_10088.txt	13210203	8f73d67a2f3b65ce237c975950be9ba6	doi:10.57745/BWU0FG	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720827
720620	Serie_10089.txt	13184166	8707f30f93aff0580826c2ec2df2c721	doi:10.57745/BIDYFB	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720620
721216	Serie_10090.txt	13163958	add72d9fc3b0c1afc886561a3a49e188	doi:10.57745/KXZG09	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721216
721153	Serie_10091.txt	13191745	821bce7ce97eff39a25e0829478340e9	doi:10.57745/UZFUXD	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721153
721135	Serie_10092.txt	13200981	082079545dcc278ce06946e54d978900	doi:10.57745/WAMG82	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721135
720836	Serie_10093.txt	13197766	228ec1ed00db0e858701c09a83a5c239	doi:10.57745/CU4RY5	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720836
721163	Serie_10094.txt	13224102	468e63eedb9c8fffd08ef4549e84d4eb	doi:10.57745/DIZZFJ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721163
721133	Serie_10095.txt	13217000	30e948ccf7e362f49db24bd641ffaa77	doi:10.57745/ZASVHA	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721133
720670	Serie_10096.txt	13223995	386ce21454ba23de30fc63ab8d9ee055	doi:10.57745/0WQNAK	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720670
721175	Serie_10097.txt	13208863	33ff6880e235858dab58970aa5c85e1c	doi:10.57745/CYSJMG	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721175
721382	Serie_10098.txt	13192155	f6935b933ec191867de2a2acf29fd786	doi:10.57745/11H113	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721382
721357	Serie_10099.txt	13179948	2f2c4292d70d8200a0ac31ff879c699d	doi:10.57745/3YU6JI	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721357
721003	Serie_10100.txt	13163783	65e05c2b592930e5f55b830b332a64a3	doi:10.57745/EIVTQM	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721003
721036	Serie_10101.txt	13175316	73281db8913a45a90537f2cc25e6d333	doi:10.57745/XSG8TX	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721036
720838	Serie_10102.txt	13209747	eec77d44c2f35cfadf04d7385baff2ae	doi:10.57745/X8BH5R	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720838
721103	Serie_10103.txt	13190879	52537f40ea21c2df3cf0f611b8760bb5	doi:10.57745/JGKQ0S	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721103
720507	Serie_10104.txt	13204654	22e3fbe74c0f34d2114db5284a09f411	doi:10.57745/RSYLZE	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720507
721034	Serie_10105.txt	13207673	ecd28d5e0dedd020740ff814c3db30f9	doi:10.57745/CMK3FN	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721034
721214	Serie_10106.txt	13198031	adbb7e6ab470eb589b8faf5dbc888ce5	doi:10.57745/HMD8JT	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721214
720795	Serie_10107.txt	13197758	ccda9f1c595600a7558f6bed29de35a7	doi:10.57745/Z7OWL7	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720795
721122	Serie_10108.txt	13207043	749bec1570008e2cc72293e4943c2e85	doi:10.57745/8Q7SEW	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721122
720875	Serie_10109.txt	13215702	79df867b8b62cfe9bd7f27872049799d	doi:10.57745/XLCOA4	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720875
721284	Serie_10110.txt	13218646	09a690c929bb3b6992ca37230ed80357	doi:10.57745/WNHP7Z	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721284
720991	Serie_10111.txt	13201096	f60801d508e9f25aebe654bf5e5c8e13	doi:10.57745/EKPLTF	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720991
720904	Serie_10112.txt	13202776	e3e4f5aefb829893cd9945f5618bf2a8	doi:10.57745/WAGJYB	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720904
720850	Serie_10113.txt	13207990	51911d865f43e8ee95a97a459208e0fd	doi:10.57745/VQQOGP	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720850
720719	Serie_10114.txt	13228821	5fa01946fec1c09a2f83e1bd55d4a4cd	doi:10.57745/UT85UI	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720719
720675	Serie_10115.txt	13197957	f5b112e2e2033e48fd9a42a5f0eced84	doi:10.57745/A8PV79	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720675
721338	Serie_10116.txt	13196082	d0ca0140eabad45df63647aee374ad0b	doi:10.57745/RX7UJ7	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721338
721258	Serie_10117.txt	13206610	14b67e478e71780c3ee19f92bcf40cb9	doi:10.57745/9AEGT4	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721258
720445	Serie_10118.txt	13207559	852f2a0fc05640c5d9cdfe372839fc79	doi:10.57745/USZHB3	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720445
720475	Serie_10119.txt	13195434	b464ff9fa270c2c11f7634d98e60917a	doi:10.57745/RJJZGF	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720475
720715	Serie_10120.txt	13204072	1d9208e270e4d995eb49aa0d07a504a1	doi:10.57745/L2BEN1	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720715
720532	Serie_10121.txt	13205335	27717d9d1797e26a444abd552a91a051	doi:10.57745/LGDHHZ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720532
720890	Serie_10122.txt	13198679	320a70d7326d5c192701c32154480e5e	doi:10.57745/EQ9CHR	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720890
720782	Serie_10123.txt	13215759	825527a986a201e50a35e4c6bb767f88	doi:10.57745/ZC56YA	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720782
720416	Serie_10124.txt	13219843	3bcd6ac17155f0464d8354feb3d3d9b8	doi:10.57745/68FRRI	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720416
720639	Serie_10125.txt	13210371	1981257da980fa82b56f16d27afe24a1	doi:10.57745/NZAWRX	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720639
721289	Serie_10126.txt	13190862	c5c6ae203fc767a460c3720237372fa2	doi:10.57745/MOZPDL	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721289
721142	Serie_10127.txt	13188372	57a48f11cedcf4e8dc52eb4ac1e66d7f	doi:10.57745/3ZNMNZ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721142
720617	Serie_10128.txt	13206206	49d41b6365a388eaf18e1ef0e8cab3cf	doi:10.57745/WWZOPE	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720617
721262	Serie_10129.txt	13193223	0916f44ab47cc6e11190753348a58030	doi:10.57745/PUVYPF	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721262
720772	Serie_10130.txt	13221310	b8a9a5bad7c061f4144515eb3c93a522	doi:10.57745/TLRWA0	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720772
721099	Serie_10131.txt	13229512	82d01fe829bc7de29ac96c3f4d76a394	doi:10.57745/0UJLCR	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721099
721253	Serie_10132.txt	13206446	7e5c4160c9e34a8770a891cc90cd16f3	doi:10.57745/QIQDEA	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721253
720558	Serie_10133.txt	13176325	a4b2d7c4ef2d2867ad67b1d9775b964d	doi:10.57745/XAPZYR	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720558
720453	Serie_10134.txt	13212650	bab14201d2c2c3f2646f6f3bcead97ff	doi:10.57745/TC12GU	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720453
721185	Serie_10135.txt	13192498	25607530bc0f91f5a3c3b1b8a87cef66	doi:10.57745/OTLHX7	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721185
720871	Serie_10136.txt	13219166	857808c2bd38e38df1ab9bb86c42a309	doi:10.57745/44YBQR	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720871
720407	Serie_10137.txt	13211597	5283a94bd922573e6181701d513e1835	doi:10.57745/6V1MQB	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720407
720570	Serie_10138.txt	13211171	3ff0d22ef43ee52f686090a38fcdf25b	doi:10.57745/KYST1X	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720570
721102	Serie_10139.txt	13221470	71cfedbd240c3bb89719abb4247e2c96	doi:10.57745/A3BG2N	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721102
721112	Serie_10140.txt	13217380	b7566b0a1588d0979f5a9535aa995e26	doi:10.57745/P7EE4L	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721112
720472	Serie_10141.txt	13191040	6fca7b01b187f0d12cee5cb53e4f8d53	doi:10.57745/JCPB3V	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720472
720459	Serie_10142.txt	13198676	eb7528d7a9dd29ca2952f1bc4e8fb60b	doi:10.57745/J7MJYE	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720459
721367	Serie_10143.txt	13222700	41e1e024a3265a56816a00a721815f65	doi:10.57745/B1ABMZ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721367
720971	Serie_10144.txt	13203820	934631b5953a3b43e36dec22a90677b6	doi:10.57745/RXFLIR	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720971
720962	Serie_10145.txt	13217699	95588c1159b9f7947a1f71ec8d0c1792	doi:10.57745/GLLRCW	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720962
721104	Serie_10146.txt	13215684	d12a4f52a76d0af61944c366c03ed51d	doi:10.57745/YELXQS	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721104
721254	Serie_10147.txt	13197618	caeb419c2184e75c28600aaeaac66197	doi:10.57745/MXIRXY	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721254
720432	Serie_10148.txt	13179833	39450cb0c5093448300ccd939c5dbc85	doi:10.57745/TMIL7U	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720432
721141	Serie_10149.txt	13177550	497b797ecdccc653c471cb37461610f3	doi:10.57745/JV1TWY	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721141
721157	Serie_10150.txt	13193118	7940675bca105e45ed1254a9be7a6981	doi:10.57745/NAEBEZ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721157
721001	Serie_10151.txt	13206613	8b9695be774e6d1831a7e6b88e110ef2	doi:10.57745/VJKC0X	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721001
721089	Serie_10152.txt	13194013	4b71efba79a8fb3a92bebb99ddbe212a	doi:10.57745/TDODDU	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721089
720885	Serie_10153.txt	13225055	e97fc59c54323e3035147d808c6900f3	doi:10.57745/TVCCYU	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720885
721385	Serie_10154.txt	13225834	b9e7e1d282c06fc8c7ce296608f72d14	doi:10.57745/12PABR	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721385
720954	Serie_10155.txt	13221524	6fe39ab3df3eba98ab81d2fabfa537fa	doi:10.57745/R6BGCF	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720954
720654	Serie_10156.txt	13192766	325d00ed4a949da72c5f3fb799619ca7	doi:10.57745/EMGV23	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720654
721272	Serie_10157.txt	13209410	6aa464e009a5c11b56fddfb1929ce7ac	doi:10.57745/LE9EJY	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721272
721013	Serie_10158.txt	13193951	9aa7dae2dec8d45228577d83c0783595	doi:10.57745/BS3OIZ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721013
720813	Serie_10159.txt	13179095	69f89bb2672a9889f8cbb43cdbdf4db5	doi:10.57745/A0WZ6W	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720813
720717	Serie_10160.txt	13187619	5d6754f51c144a6a8c83b9088bf73b7c	doi:10.57745/WRRX3U	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720717
721064	Serie_10161.txt	13213071	ec0ff64a7e54f8e1ef8ff2d29dcc528c	doi:10.57745/9KCPV9	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721064
721072	Serie_10162.txt	13201521	eaee3045dae1460e4754955eaaecced7	doi:10.57745/UQ4XGH	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721072
721092	Serie_10163.txt	13217692	473bc2e1b6ec505a144ce87baae70839	doi:10.57745/MVSLFL	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721092
720746	Serie_10164.txt	13206032	0cd3037cadee2dd17ccf685673cf873e	doi:10.57745/FQ9DFE	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720746
721000	Serie_10165.txt	13203603	e97b9b063c93f26786dba30d2ab386fd	doi:10.57745/H3WXFK	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721000
720556	Serie_10166.txt	13189362	807451bbc47ec987f21ffb8089720335	doi:10.57745/GB6KGP	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720556
720583	Serie_10167.txt	13186842	07a28c1f1897861c88d88f1839196ae8	doi:10.57745/SSFPPE	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720583
721106	Serie_10168.txt	13207896	9cc6621224984fb199d7c140748d74d8	doi:10.57745/PDN6SP	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721106
720801	Serie_10169.txt	13228919	8873c9d3a5ac3aa8cccc5c77fe3cf9ba	doi:10.57745/U8EJYF	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720801
721197	Serie_10170.txt	13197696	e6100f7e1d3d7ad7a809abd973b6a22a	doi:10.57745/F5BHL1	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721197
720664	Serie_10171.txt	13217307	fb472c12307ec8add7b5d52566671701	doi:10.57745/O0WANR	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720664
721341	Serie_10172.txt	13191383	13b8cf85d171fc504eec82cd983e27a8	doi:10.57745/XOPI84	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721341
720567	Serie_10173.txt	13216176	18849a2f16e62083563ab6b03e4df505	doi:10.57745/BOIPR2	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720567
720996	Serie_10174.txt	13195081	50c87bff5f1a5fda1a251a70619e47a2	doi:10.57745/GRWC7K	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720996
720536	Serie_10175.txt	13183131	d7215932670f5a7f4c22af6119979216	doi:10.57745/VLHCOC	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720536
720825	Serie_10176.txt	13197766	e7888aec20cbaf959ba4009bb90a5228	doi:10.57745/KDMXU7	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720825
720577	Serie_10177.txt	13202858	71f2ce01276215fc588874641188581d	doi:10.57745/ZR5SDK	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720577
720553	Serie_10178.txt	13187466	5423e4623e7a30a788e2c1380a983650	doi:10.57745/LJHOMH	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720553
721074	Serie_10179.txt	13204741	12bc7cae41c8c694a2abec34e7845cc8	doi:10.57745/RNC2YJ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721074
721091	Serie_10180.txt	13205966	83489dd983ff79ffabfc1fc49d760124	doi:10.57745/2BGCYC	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721091
721165	Serie_10181.txt	13223119	66fb1002600c251e83e77b4b1dde8e69	doi:10.57745/ITIQ8Z	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721165
721387	Serie_10182.txt	13199171	d692761b3ce1f06bd321eee323b7dd5d	doi:10.57745/YJ9ZFE	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721387
721348	Serie_10183.txt	13170268	66b92c26f2b1b94117b4a48d6de5a020	doi:10.57745/3KANS0	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721348
721363	Serie_10184.txt	13181901	5742842ef86a98af1b313192f73b7b52	doi:10.57745/HPQYD5	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721363
720576	Serie_10185.txt	13208575	186e519044b9de148be07ed2036056e8	doi:10.57745/YRJRB0	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720576
720423	Serie_10186.txt	13209297	e5771c9434c1aa4d20d5b2ee5c8779d9	doi:10.57745/L3WAYN	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720423
720903	Serie_10187.txt	13194016	df33d77f3e18dc5337840bb237c9a541	doi:10.57745/HB3MOL	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720903
721301	Serie_10188.txt	13203382	97e5813472de9d6ed992c02acf9e1f2e	doi:10.57745/DTF0LX	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721301
720519	Serie_10189.txt	13207668	a96abce9deab6f49e18de8cb11111de5	doi:10.57745/9MRH6M	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720519
721120	Serie_10190.txt	13188502	bfb6575a2b4c818d9325cb77888bcacc	doi:10.57745/LBUJQB	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721120
720530	Serie_10191.txt	13196315	36c40f8532a544d12c6ce1a21327a0a1	doi:10.57745/TJZ4ZZ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720530
720953	Serie_10192.txt	13202187	45bc94a654acda605f1228445ae08f40	doi:10.57745/KZHWLO	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720953
721127	Serie_10193.txt	13213913	a23c49c1c9c4546f918df676d73c0c72	doi:10.57745/B1NE07	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721127
720495	Serie_10194.txt	13197964	c92e7af2fde8b502ac5bbffd44e24456	doi:10.57745/JCHOC6	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720495
720581	Serie_10195.txt	13188300	3657c57343677394a67b28a77d91c5f7	doi:10.57745/L9IVP6	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720581
720865	Serie_10196.txt	13183568	591a570ee09679e1808c597254e67dfa	doi:10.57745/QYUDWO	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720865
720913	Serie_10197.txt	13204404	f9d8000ac3256416df5ccccf67dd2004	doi:10.57745/ZBIBLK	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720913
720780	Serie_10198.txt	13181033	896136e872074fcbb9c90349868c00a6	doi:10.57745/HQJLZW	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720780
720673	Serie_10199.txt	13212124	016c9a615a8b17364c1c593cb0dec6eb	doi:10.57745/G113CF	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720673
720595	Serie_10200.txt	13181222	6813441246240c3f6928dd298330624f	doi:10.57745/U12MDB	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720595
721303	Serie_10201.txt	13171762	7866b2bec55bd8a43be7d5f2a1024ad1	doi:10.57745/O9XMJN	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721303
720622	Serie_10202.txt	13221691	5da6f68e97f8d67e2d7d1784a5bc6c7a	doi:10.57745/IPRS6V	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720622
721235	Serie_10203.txt	13218064	db25b77a9c126944e79fc2d8407a6b58	doi:10.57745/BNQ0XX	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721235
720481	Serie_10204.txt	13219596	0dbecf3d3a0149db9a5fe9c0996618ed	doi:10.57745/BGBSED	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720481
720516	Serie_10205.txt	13196897	877e66e05faf97fda862a45e9b72d392	doi:10.57745/ODWR7Q	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720516
720494	Serie_10206.txt	13201824	8f6f821fae366480565c8e8664105d11	doi:10.57745/Y2V368	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720494
721377	Serie_10207.txt	13182671	1916ec84b9e6a0d3318272df3eb0d30c	doi:10.57745/YDCLLY	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721377
721115	Serie_10208.txt	13205073	a53826eab10795a011d7dc6c65ae33c3	doi:10.57745/TTYGRM	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721115
721334	Serie_10209.txt	13223242	7419fc63b69d0145ea5b9b15e77234c3	doi:10.57745/RGY5DZ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721334
720816	Serie_10210.txt	13210134	22bdd4eb0e53781db855e54219b39156	doi:10.57745/HHGBBA	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720816
720881	Serie_10211.txt	13187023	3e0f9b1e50fa8da9199b4141bcaa3203	doi:10.57745/FYZLM0	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720881
721131	Serie_10212.txt	13210433	8a59c0f33415b0a741338c67952ae9b6	doi:10.57745/HXVELQ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721131
721048	Serie_10213.txt	13218534	faa3fac3503497ae3d624902de3c44e2	doi:10.57745/7MQ1YY	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721048
721125	Serie_10214.txt	13218030	65310950aae579a5ef60ec498a6246ad	doi:10.57745/UOTLPI	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721125
720552	Serie_10215.txt	13197187	2e00880dc4aba7533a00decf543f99da	doi:10.57745/X0JP0T	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720552
721236	Serie_10216.txt	13204019	9d93d6e5a77c295bfb4ef7ecd9ab07a4	doi:10.57745/23PPQM	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721236
720783	Serie_10217.txt	13211181	75bcc2acf8741f3c9a00a7da7895fbe4	doi:10.57745/LI7TIE	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720783
720800	Serie_10218.txt	13190560	d2e20b6a0427e9a8d506563ca131b5b3	doi:10.57745/IM6LVX	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720800
720511	Serie_10219.txt	13177506	c5a59ad046c522c2d14c3fce8fde0cec	doi:10.57745/89L86Y	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720511
720436	Serie_10220.txt	13193661	3f3bab28172530ea2cf4ed7766a7e0f9	doi:10.57745/E8L0E6	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720436
720776	Serie_10221.txt	13176176	06f17dedbee9982095be432be6b44e1e	doi:10.57745/RZWBQS	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720776
720643	Serie_10222.txt	13182874	1bbaa0da62f6492d1726cad6b1589a97	doi:10.57745/RCYS4X	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720643
721323	Serie_10223.txt	13212923	4e30583e37afe87afb555acca24b2b5b	doi:10.57745/4DVSNF	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721323
720447	Serie_10224.txt	13205761	7891a9e74f246f332b55c22b5c494288	doi:10.57745/PFSWAS	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720447
720909	Serie_10225.txt	13209052	be1abb0753c12fa63b22c9f0de5fd9ec	doi:10.57745/TAKKMH	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720909
720535	Serie_10226.txt	13202804	781c22cf9a09f7bdc1668fb2c836f298	doi:10.57745/GESOHG	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720535
721324	Serie_10227.txt	13196904	5fbdb1b3213d1cd0d719409d2536abd2	doi:10.57745/WENFZJ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721324
720848	Serie_10228.txt	13210543	f68029495c5dcba19f0cf7078d465506	doi:10.57745/RGDF4L	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720848
721015	Serie_10229.txt	13206792	8c5c334051a844a1a196570aa95ca815	doi:10.57745/XHULJM	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721015
721116	Serie_10230.txt	13222974	c1ae4004d8fd1d6d31db17ee6442287b	doi:10.57745/DRWCFW	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721116
721347	Serie_10231.txt	13224914	9c995f7969ee734891bb03c4ac883439	doi:10.57745/2F1BZV	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721347
721364	Serie_10232.txt	13209941	d3fa1347f6e9238339682cdd6eb162e8	doi:10.57745/QGUAE1	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721364
720901	Serie_10233.txt	13238296	ea19c78f93fdbe8e23290ef0e1e7ecc1	doi:10.57745/HFMGGE	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720901
721168	Serie_10234.txt	13193754	955f6657d69d116613d3c5c029adb052	doi:10.57745/XCVV0Z	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721168
721011	Serie_10235.txt	13185256	e0c62a2c588b0aeb8cc07d1a8350541f	doi:10.57745/QGA2IK	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721011
720655	Serie_10236.txt	13188515	cc0e1c52492c556d37db75f4866619fd	doi:10.57745/5OHBTP	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720655
720732	Serie_10237.txt	13203477	5e6056167a59ce0f6303acb499cb1cfe	doi:10.57745/O3MPDD	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720732
721073	Serie_10238.txt	13226046	3a44c29648b083e7490a9ca0cd4a5739	doi:10.57745/NHQQNH	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721073
720756	Serie_10239.txt	13200442	318c70f2d8cf6dd5a939bc17e5786029	doi:10.57745/ZUVXYI	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720756
720785	Serie_10240.txt	13215067	fab2a2d6664099fb4651b3c3a2bec621	doi:10.57745/V1OYTV	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720785
720886	Serie_10241.txt	13203263	25dcba45da78f9ec23489f353bfd0a5e	doi:10.57745/ROKSTD	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720886
720766	Serie_10242.txt	13206099	cac146aec82b9a83e1cb8e94b306b1e4	doi:10.57745/GTATKB	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720766
720454	Serie_10243.txt	13202144	19aa65c165247320461ae692912f94f0	doi:10.57745/JJIW1H	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720454
721039	Serie_10244.txt	13207074	774707ab324b29ddc57c8db98f39876f	doi:10.57745/XTUMWE	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721039
721124	Serie_10245.txt	13193259	1971623e5592f1cc83d5700f32852b8b	doi:10.57745/NG1FX1	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721124
720740	Serie_10246.txt	13206217	7d896c856239aeb61754e161619e2e27	doi:10.57745/OJFSLB	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720740
720888	Serie_10247.txt	13228687	569999694fcbbd46666ac7d589ffed5f	doi:10.57745/F2UNMK	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720888
720460	Serie_10248.txt	13207655	9ac7ab13ac7a22d1722c29f44c99760d	doi:10.57745/ZKSA50	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720460
721088	Serie_10249.txt	13202876	20daf50a79fc61b8dc1306c5167ea907	doi:10.57745/EQUM3A	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721088
720627	Serie_10250.txt	13199083	504fe56bbfbebd685dcd83ae30096ce2	doi:10.57745/TYVBQY	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720627
721083	Serie_10251.txt	13194228	c4a5131e43ddcbd2f6ae45c0cb480de0	doi:10.57745/GTJ9RW	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721083
720505	Serie_10252.txt	13182449	48c7856089a569312b1e32edecb3a8f1	doi:10.57745/X8MFJN	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720505
721281	Serie_10253.txt	13223132	a4016f1c824ac2017f17fa6318f886ab	doi:10.57745/YHSSNI	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721281
720993	Serie_10254.txt	13187919	93ada063ca5f4b174c0ddab5ef157d54	doi:10.57745/KVABKV	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720993
721213	Serie_10255.txt	13191822	c0756b3560e7eafd3041c855b75fe0bd	doi:10.57745/FQ84XV	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721213
720960	Serie_10256.txt	13186363	3d514d4d7edad4895a433a7442e6d7de	doi:10.57745/XLDOLH	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720960
721221	Serie_10257.txt	13189411	ea0ec2d83d33feacbce13d3a7794b726	doi:10.57745/IIMIUB	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721221
720660	Serie_10258.txt	13193801	a04f38982f89e4d19c0b29a080b0e413	doi:10.57745/9CMKJZ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720660
720404	Serie_10259.txt	13187584	95b4691be76c11366d6f23b0cd8ebd13	doi:10.57745/MNZWP2	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720404
720784	Serie_10260.txt	13240610	3d56df36761cb8dd51ea3bc28e169f62	doi:10.57745/D4OIKP	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720784
720589	Serie_10261.txt	13226833	ae8c8d9ac426e1fa64db9c3660a371e4	doi:10.57745/4ASG5U	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720589
720616	Serie_10262.txt	13214279	770b615cfbdbe92cd3fc2c82c529bd0a	doi:10.57745/I2CYI4	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720616
721276	Serie_10263.txt	13207721	96f673e79d7eaf62458456b61c3fdff3	doi:10.57745/S9MVYS	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721276
720403	Serie_10264.txt	13191278	9dc53decb55a31bcd47675200a99ed09	doi:10.57745/VGFCXZ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720403
720431	Serie_10265.txt	13215720	c4d9e16d931cd51ce7de6a548a22f2f0	doi:10.57745/NC22B3	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720431
720523	Serie_10266.txt	13214715	82fa8aa81302189753c0ede3a48add67	doi:10.57745/FYIGGA	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720523
720727	Serie_10267.txt	13203280	42cbfdf9e656bbec1f7cdcf7c52019a9	doi:10.57745/1UU6YD	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720727
720579	Serie_10268.txt	13180088	b1edf1e05b4e90a8d0632e42939ff38b	doi:10.57745/13SXVC	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720579
720990	Serie_10269.txt	13207119	314f3be70b254136c0e7a289e7ea49c6	doi:10.57745/KE7EVO	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720990
720751	Serie_10270.txt	13229799	6cecf23568e52d49ccff35f3ceaba23e	doi:10.57745/DV7FRD	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720751
720450	Serie_10271.txt	13182537	c4f3667e951493fe67539a76d5c93458	doi:10.57745/17KKTI	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720450
721369	Serie_10272.txt	13200522	f1002f825c654e145da309164551374d	doi:10.57745/F61CHJ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721369
720662	Serie_10273.txt	13211294	2aed8ebf977bc521391927a07457a0aa	doi:10.57745/XA2KSJ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720662
721307	Serie_10274.txt	13188928	b34ba01a18fdc6bd4e30877b292b987d	doi:10.57745/1ZTWU6	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721307
720433	Serie_10275.txt	13181288	760a5b41f78aa798ecc4121871c2fae1	doi:10.57745/DWZTJT	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720433
720694	Serie_10276.txt	13203508	7b8d324037ea389de0421a6abb8cda4e	doi:10.57745/YHTE3Y	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720694
720919	Serie_10277.txt	13195746	66af55115a1901991e4bb92f94f96f04	doi:10.57745/QVRTKZ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720919
721005	Serie_10278.txt	13205341	59148658e171538b24598bf41f8ac2c8	doi:10.57745/UREVWQ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721005
720731	Serie_10279.txt	13203209	324e40d0e0545b0a441fe3bd5b938b24	doi:10.57745/XAMEMI	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720731
720496	Serie_10280.txt	13220433	648b15765358682a0f494e22bd653472	doi:10.57745/JGDMFM	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720496
720529	Serie_10281.txt	13202961	89fc9edb4e8421cac625153eb6bc496d	doi:10.57745/VSFXVB	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720529
720765	Serie_10282.txt	13193787	38aae9a3ed4fc6ad92ab012e95217da5	doi:10.57745/SGXHE6	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720765
721267	Serie_10283.txt	13193697	508473a9035be3db19ca41a0d66a1489	doi:10.57745/NJBABR	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721267
721285	Serie_10284.txt	13206158	25763a4c129699d27adfd9127c94da06	doi:10.57745/B66DC7	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721285
720743	Serie_10285.txt	13204285	4faf050e9dc6968fbdda9e1331a6fc08	doi:10.57745/PTY7FN	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720743
720709	Serie_10286.txt	13211248	85b8d9bc3c066d5395ab9ecb45a55aba	doi:10.57745/KWVTWR	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720709
720659	Serie_10287.txt	13210545	cbb26ef060a4ce7cbf84d628f7bb5e37	doi:10.57745/JTJAQY	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720659
721033	Serie_10288.txt	13177544	1d33d95d68efdf4eda22bdc04c1f57b7	doi:10.57745/C3MNP5	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721033
720462	Serie_10289.txt	13186735	2407f816cbc1db4fbb1307d4e7d500f9	doi:10.57745/LJDO5V	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720462
720973	Serie_10290.txt	13180961	46334292384c65bb1e15275a990697ed	doi:10.57745/CKY7DU	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720973
720674	Serie_10291.txt	13190832	020094ec7173c7ade4294a61c7854273	doi:10.57745/JKOTZV	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720674
721376	Serie_10292.txt	13194953	ab92b9e072f18d5d55bb31c2b0026b3b	doi:10.57745/UFTOFJ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721376
720964	Serie_10293.txt	13213040	4ac11a69133766b9203d92cf58bb2004	doi:10.57745/CCUGOS	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720964
721147	Serie_10294.txt	13186080	3d1215259902917acaf4c6af953357a0	doi:10.57745/WNPQ8I	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721147
720585	Serie_10295.txt	13225845	c1784571a33e8ad59a97b76ce87bad62	doi:10.57745/3PZXF4	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720585
720582	Serie_10296.txt	13194352	ed8bf973a8846c04806c7f9aaab58858	doi:10.57745/62T1AB	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720582
720917	Serie_10297.txt	13204243	21f5ff22025e9d7cbe177791fac20d7c	doi:10.57745/AIJZ8J	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720917
720599	Serie_10298.txt	13197607	c32f86e7ffbc3cc9974066ba168de626	doi:10.57745/8DXRJC	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720599
721136	Serie_10299.txt	13216300	e0b02d8fb89c4f782e73b1bc98034ade	doi:10.57745/GSBALY	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721136
720586	Serie_10300.txt	13215853	840db3e43a4ab9a78007a6df4256d9ac	doi:10.57745/R77XIE	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720586
721340	Serie_10301.txt	13208073	8267ba00578091171607f731dc622d75	doi:10.57745/RE0VKY	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721340
721215	Serie_10302.txt	13226150	3622ce3e210625b11a7c08bf54847f19	doi:10.57745/MYOF7T	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721215
720793	Serie_10303.txt	13209530	c5a78c51abdba67e9a30d7be8bcfbb33	doi:10.57745/DFTUBV	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720793
720613	Serie_10304.txt	13198648	f83ebb08f3918eb78210f51f8460792e	doi:10.57745/UQMW45	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720613
721222	Serie_10305.txt	13210667	6ca216889d41a657ef67ccdd27b93af9	doi:10.57745/LUXP73	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721222
720506	Serie_10306.txt	13177680	b4c258972983c315bbb52aa828701661	doi:10.57745/HGGKTJ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720506
721339	Serie_10307.txt	13189799	50b4e8989419cc644db7ac8b190a891f	doi:10.57745/JOSQKY	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721339
721032	Serie_10308.txt	13218024	12dd5ef2c3324eb31939e0c6f1608dea	doi:10.57745/RWX4NW	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721032
721191	Serie_10309.txt	13204784	b2d201ded24b2a0bb233d04b62c2dd60	doi:10.57745/IREY1M	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721191
721278	Serie_10310.txt	13226238	face3c81f731b00b299966b202479f4c	doi:10.57745/SCNQNB	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721278
721014	Serie_10311.txt	13185177	538c92314b3e9cb23d453a6f9ef10bda	doi:10.57745/ONFFRL	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721014
720665	Serie_10312.txt	13176090	627719b17b68380b51bf0f32c44605d2	doi:10.57745/LLWEYX	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720665
720685	Serie_10313.txt	13198108	58211105277e94861e6965f53be36df4	doi:10.57745/EF7TYZ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720685
721368	Serie_10314.txt	13207616	edc65842cc02ff09b28cbc68bfc6219e	doi:10.57745/Z8YD6C	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721368
720656	Serie_10315.txt	13236010	bad373c7487cb35af11fdf9924e0707f	doi:10.57745/9BDIWE	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720656
720490	Serie_10316.txt	13210625	0690f81826cee1049963397af5dcaa86	doi:10.57745/PIR40U	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720490
721379	Serie_10317.txt	13224617	b1b916507eecd3c93d328a1d2efd6ec9	doi:10.57745/CAAKMM	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721379
721343	Serie_10318.txt	13205001	21c46ea5a714b2e867c0768f821d1ed0	doi:10.57745/GABIOX	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721343
721055	Serie_10319.txt	13183991	f44085fae3388eee21aaa58c895b8524	doi:10.57745/NUI7MP	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721055
720729	Serie_10320.txt	13201231	cb7d7a2f7b7744ceb403569e395b94b8	doi:10.57745/LBYAWO	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720729
721019	Serie_10321.txt	13199509	b0c16d30da18db9346e9d05945db7068	doi:10.57745/SXQMDA	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721019
721212	Serie_10322.txt	13195877	9407b7c3ab9905e44aa5b6a61827d5d9	doi:10.57745/OSMIYK	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721212
721362	Serie_10323.txt	13187313	b670a891fdde23a902868f5ab3cef38b	doi:10.57745/6TVWG5	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721362
721044	Serie_10324.txt	13213692	410c9d80ee1f9a1cdc1bf84c689424b1	doi:10.57745/MUGWLA	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721044
721358	Serie_10325.txt	13199062	2163f765b206e9f8f89a75ab4d139cb1	doi:10.57745/NHBEVO	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721358
721041	Serie_10326.txt	13221956	7a25c02f1f9ae5772e1d3ae29fb5f520	doi:10.57745/YRN8PR	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721041
720747	Serie_10327.txt	13208290	49669ea54475c053d53d0730379283d8	doi:10.57745/RXYCXP	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720747
720724	Serie_10328.txt	13203885	74e44d2bb16cbd7c3ddd6d2635fae965	doi:10.57745/HXMTEZ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720724
720737	Serie_10329.txt	13184323	5d4c5c0ffd98f3e7c7bc2e792c7ba2d7	doi:10.57745/ZABLTH	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720737
720720	Serie_10330.txt	13177210	388929e427f81a06969ff247e86b9675	doi:10.57745/LMKVWX	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720720
720837	Serie_10331.txt	13177088	c5527ecd13b36808fa2b3c65d859f99a	doi:10.57745/JZM4IV	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720837
721161	Serie_10332.txt	13167568	33d61a4fce330e8cd7804696c97f5908	doi:10.57745/6LMEHB	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721161
721263	Serie_10333.txt	13208452	5c81cc1c3466261cf1c8abb88e2bbde9	doi:10.57745/2MXX99	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721263
721060	Serie_10334.txt	13228305	80db666fcabf9d25e3c198be6fb99100	doi:10.57745/QL2Y63	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721060
721186	Serie_10335.txt	13209496	cd8ef00d4f05241ca7aa691d211fa6ae	doi:10.57745/Q6QRE3	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721186
721381	Serie_10336.txt	13199348	c598d4898a5cd6082018ce5b7dbb5acb	doi:10.57745/LMFIC3	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721381
720597	Serie_10337.txt	13183285	5ec2ffb8d41e01b4c3b3bb60a991c481	doi:10.57745/P4OUSC	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720597
721204	Serie_10338.txt	13200058	dd47d2c1a1a43009909d1c4a5b4c494d	doi:10.57745/PFTU0C	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721204
720415	Serie_10339.txt	13206501	046d56e721f79de2f7f64cebf42d46fa	doi:10.57745/TDZV0P	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720415
720722	Serie_10340.txt	13203958	b108f1174cbd2ca418e0ebc39eb6cc21	doi:10.57745/WA9VML	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720722
720648	Serie_10341.txt	13228955	9440eb419b2251dcefb27c501529ddb5	doi:10.57745/BVD9L0	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720648
720896	Serie_10342.txt	13194996	10a67a764457ee0bef07d3af29a91a72	doi:10.57745/PEJIBL	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720896
721355	Serie_10343.txt	13207054	1264c6171b9b6889f8b3d851875fdbcc	doi:10.57745/AYGQEO	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721355
720995	Serie_10344.txt	13204085	c7a0bedebb9169df34dcad9c0791ba84	doi:10.57745/KDMWSZ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720995
720686	Serie_10345.txt	13198369	e72aa66498f421403068ff1d4b751b0e	doi:10.57745/08DLTY	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720686
721240	Serie_10346.txt	13190083	45b72b550e5d516eab3c1ae65db74a54	doi:10.57745/WODRDJ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721240
720452	Serie_10347.txt	13183383	741074757d7fa6e3bcbb7fb4a87d2b23	doi:10.57745/RQ6UWQ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720452
720695	Serie_10348.txt	13195718	fa2fed327cf5f7a9f29c77235d80f2af	doi:10.57745/ERBPVU	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720695
720533	Serie_10349.txt	13245568	6568dac480599f14f62a50c03d2ddc1c	doi:10.57745/OZMROM	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720533
720653	Serie_10350.txt	13207788	4e602afb554c01a63da50c13c72e933f	doi:10.57745/U08TWT	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720653
720607	Serie_10351.txt	13183413	9dbefec1bb7978d44e42f019890caeca	doi:10.57745/3OWCMG	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720607
720930	Serie_10352.txt	13204530	bffaac506e2ab5b38262789655330221	doi:10.57745/IUGKEC	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720930
720986	Serie_10353.txt	13212261	a19439e7299d33df5bf4f96fe3899ef4	doi:10.57745/XLOCZO	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720986
720687	Serie_10354.txt	13218847	f2e047dc23c4f97ba56941009de239ef	doi:10.57745/T6LHZB	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720687
720565	Serie_10355.txt	13184758	6bb9c882c8cc18951c0ced28625d5433	doi:10.57745/PWJNOK	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720565
720414	Serie_10356.txt	13169998	f514b12ecabcfee8f72d353428529dc5	doi:10.57745/CPQXGW	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720414
720625	Serie_10357.txt	13227644	27b1fd5ecbcf2deefce87d3dcc556249	doi:10.57745/LT2PJA	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720625
721069	Serie_10358.txt	13216253	58ee141f2f962e1b0cbfbdefd898dc16	doi:10.57745/JUIMLM	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721069
721126	Serie_10359.txt	13204548	b287f51955fc32a20682512c9fac02f3	doi:10.57745/M312CN	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721126
721344	Serie_10360.txt	13220236	d4a7e4bdaec35760af9c0c459d92a21b	doi:10.57745/IIXBBG	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721344
720534	Serie_10361.txt	13189165	aa1c61d8f8ad79175cbd73239418aa15	doi:10.57745/RC22PV	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720534
720435	Serie_10362.txt	13190126	3e9d14fffb0ac6419f63a4bc9a4a9b63	doi:10.57745/ZMEYMA	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720435
720559	Serie_10363.txt	13203145	0551c743c3157f1f4c2eda75e246bdae	doi:10.57745/KOAEFO	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720559
721181	Serie_10364.txt	13189708	63f42b7fb0aa65b2e542ac22d026ede6	doi:10.57745/Z3QWQB	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721181
721075	Serie_10365.txt	13206566	2707d64e5ca447ec5b39c341da0dd312	doi:10.57745/RTWDND	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721075
720703	Serie_10366.txt	13226255	254296ec8d74bae3058b332a27a080c1	doi:10.57745/VI20CD	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720703
720957	Serie_10367.txt	13193070	660853981c27036410c9633169601484	doi:10.57745/DMPAQO	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720957
721306	Serie_10368.txt	13201887	2e2077fb4c33a2796edf14d76e969c72	doi:10.57745/GGHCFH	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721306
720438	Serie_10369.txt	13180070	afe6bb06803fbf0410b1a4bca0020a2c	doi:10.57745/L4BLIT	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720438
721077	Serie_10370.txt	13192685	c79fcd00a52c11723712a6372f372b73	doi:10.57745/PVKQLE	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721077
721095	Serie_10371.txt	13182649	965e0b3cd420c9d21cddc340897ce2de	doi:10.57745/52PKUC	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721095
721239	Serie_10372.txt	13205346	b6be3582e53631c4042efa88c7244669	doi:10.57745/0QJMWL	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721239
720498	Serie_10373.txt	13207463	1c62ab8f0ef98598fb5bd00cdfde0c27	doi:10.57745/OBPK1R	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720498
721182	Serie_10374.txt	13203679	4944a6a4097ce5074d0c95e256122378	doi:10.57745/WBACPU	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721182
720860	Serie_10375.txt	13200123	6978500442db83a5a9349c1d51b79c9d	doi:10.57745/MTVEGU	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720860
720584	Serie_10376.txt	13197786	0cac6123f0cb7200768eca141c8d49e5	doi:10.57745/KYXGVS	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720584
720634	Serie_10377.txt	13190048	58ffec25e999ebfaafe28a0786434f57	doi:10.57745/RQ86UC	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720634
720706	Serie_10378.txt	13206265	a5254c0571b78d04993672c5ba1c696b	doi:10.57745/MCU8L3	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720706
721170	Serie_10379.txt	13177470	6e75f9b05dd863671bd7acbea6f2b7de	doi:10.57745/LD1JQG	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721170
720483	Serie_10380.txt	13186316	716e0ba1968fd28ca3cd20a5265373e6	doi:10.57745/7Q7WGJ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720483
721353	Serie_10381.txt	13201275	62d0caa93b24fe48e5e4d9806b8554f4	doi:10.57745/EWPJCU	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721353
721225	Serie_10382.txt	13225692	e6c8a9d381f228eec0c7e8399ee79749	doi:10.57745/JNT4PR	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721225
721207	Serie_10383.txt	13197539	cd62c0e5decbc4326c22cc8f331a5188	doi:10.57745/TH4BVK	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721207
721030	Serie_10384.txt	13179185	11f12d87d48ef6c38e2177e78a5c46dc	doi:10.57745/MSM4MF	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721030
720841	Serie_10385.txt	13206646	08f4c4a4740c1c1dbd84ae5d0590423e	doi:10.57745/9XXGE6	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720841
720668	Serie_10386.txt	13230332	e259e073ec85c91a0375ca0bca110c3b	doi:10.57745/Y60LPI	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720668
720411	Serie_10387.txt	13214166	bf0f8403f692ec15e2d49d7854a7fbc1	doi:10.57745/WN9T8Y	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720411
720770	Serie_10388.txt	13217334	c7689656e022c0462bc051978881a9f0	doi:10.57745/P1X04Q	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720770
721173	Serie_10389.txt	13192817	74f20c4448b50e7f97aed475259da71b	doi:10.57745/YAIFPL	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721173
721292	Serie_10390.txt	13217670	076c72bb78a2f839da3b58405921b7cf	doi:10.57745/VNGEWU	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721292
720762	Serie_10391.txt	13196126	8f27aafaef4f415f21dbf64836fcb801	doi:10.57745/MR5KH3	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720762
720982	Serie_10392.txt	13218472	57c551c3c2a5fab8ef09bcf6bc48d513	doi:10.57745/VHXUSU	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720982
721066	Serie_10393.txt	13197176	2d5b7c648215da5a62978d4b4e9a3b6a	doi:10.57745/P2KDZR	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721066
720778	Serie_10394.txt	13204898	a681e9ec9f50b4672816ff5283e0ef9c	doi:10.57745/KQMQGC	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720778
721008	Serie_10395.txt	13229504	bfacb2aa23570b4e98866222b8c516d8	doi:10.57745/I6YCNY	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721008
720466	Serie_10396.txt	13216569	628aba4107b428bf13c176ed1f609609	doi:10.57745/QFQXOJ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720466
720758	Serie_10397.txt	13178827	3645bf9784e0ed1783502c21d65189c9	doi:10.57745/QPGS5H	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720758
721386	Serie_10398.txt	13183503	c4aba77738bd0a348dd57779b2bb9903	doi:10.57745/RMBJNB	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721386
720891	Serie_10399.txt	13191589	0607cdaa0e2f37e675386e1895f4dc7a	doi:10.57745/2IMEGR	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720891
720941	Serie_10400.txt	13201700	473c27e9053f8ab0f924a03d598c6cf1	doi:10.57745/EOVQOK	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720941
721243	Serie_10401.txt	13238232	e1acec77f99c7c16356aba61c4f95138	doi:10.57745/JPC93K	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721243
720777	Serie_10402.txt	13189737	09fb0b458d87bb43bf5a1d964639f173	doi:10.57745/E3IXQ4	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720777
720944	Serie_10403.txt	13200559	2b447760c18e4b2c2a8a8c2ee7db6671	doi:10.57745/KBS0WJ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720944
720807	Serie_10404.txt	13179723	7cd7688e19b7c02b9ca1827274114008	doi:10.57745/E4N1VL	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720807
721283	Serie_10405.txt	13208220	7247bf144ae6bbcf2880061f74395d15	doi:10.57745/ZTC1LB	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721283
720810	Serie_10406.txt	13219330	90e5dc4ebf7da43b284280fa2cf5170c	doi:10.57745/X4YLLJ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720810
721016	Serie_10407.txt	13208789	ed661214491a47ef14b289c7bfe936c2	doi:10.57745/VJYIZX	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721016
721294	Serie_10408.txt	13220244	bfce6807e94667a525da9f14bbb03732	doi:10.57745/UUBMUI	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721294
720676	Serie_10409.txt	13214920	647a1258fb280b69dafe5cd9f269d5f1	doi:10.57745/EER3KY	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720676
721230	Serie_10410.txt	13194625	33e5b3f40eebdf33779a14dcd49f7d6c	doi:10.57745/SMGU6Z	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721230
720994	Serie_10411.txt	13177446	930aa354fa8effcc4d03e10c39357304	doi:10.57745/RIMIBE	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720994
720626	Serie_10412.txt	13214800	93da896bfab005e5d44166f7f2cb76c6	doi:10.57745/MQ4EOD	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720626
720928	Serie_10413.txt	13186178	0d3ee2b5e6c3042e78f06b2f0009f9c5	doi:10.57745/LBB4FV	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720928
720478	Serie_10414.txt	13196213	1ae2de2ae70cdc729e4f69cbf4b429e1	doi:10.57745/UMZGND	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720478
720455	Serie_10415.txt	13207263	3f6e49e27ef29d6347dc450d92141864	doi:10.57745/KVFUF2	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720455
720855	Serie_10416.txt	13240249	620a2fd634a5ab15092f2dddb8946a87	doi:10.57745/FPA89R	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720855
720921	Serie_10417.txt	13208973	2e8d51baadbe305d3faf85526cc9901c	doi:10.57745/MOWKMM	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720921
720867	Serie_10418.txt	13206799	948e78939ef3f880994252b5ea270a0f	doi:10.57745/ALNKG3	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720867
720527	Serie_10419.txt	13242860	33d9b352df885de1c6b2ee4eafd15797	doi:10.57745/10M5GN	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720527
720832	Serie_10420.txt	13183487	f214bcfd22a4df080d82c4e231ba014e	doi:10.57745/QNODHS	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720832
721335	Serie_10421.txt	13193129	75285ac012c05f16502467f3907438c9	doi:10.57745/JZMRGS	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721335
721006	Serie_10422.txt	13193445	d4cb681f599910bc7a7526aa97f838bb	doi:10.57745/51F8B5	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721006
721274	Serie_10423.txt	13213258	e92032c3d04522bb17476038a712200f	doi:10.57745/2A3YTL	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721274
720439	Serie_10424.txt	13232249	e8337f75f47a7a830a38c39353a259fb	doi:10.57745/3S5RCA	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720439
720951	Serie_10425.txt	13215948	4e2c7b1c6d542694a6eb97b84a98c960	doi:10.57745/WA67OS	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720951
720692	Serie_10426.txt	13209566	b6a81307bd32aa54ed119aabec1be2bb	doi:10.57745/ESSIWK	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720692
721326	Serie_10427.txt	13193572	97e32aa9c8cbf5cbba48a2300ca0ec54	doi:10.57745/3IHVKN	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721326
720922	Serie_10428.txt	13195946	593098ec50e103bf014497709a1374c6	doi:10.57745/DTY1CZ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720922
720711	Serie_10429.txt	13219699	30671a54e24faf06f082108c36d03eb2	doi:10.57745/KZWRXX	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720711
721327	Serie_10430.txt	13203116	d88f70f43ed23bc832181df377e9937c	doi:10.57745/NBQYIB	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721327
720713	Serie_10431.txt	13206190	1ca9d8676a9243155c2489f0d4aa78a6	doi:10.57745/7E3SJN	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720713
720632	Serie_10432.txt	13238459	794280388f3e2f7e89ca9932dcd6d815	doi:10.57745/Q6IRMX	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720632
721309	Serie_10433.txt	13230494	2d55cf5fbdbc29ac9e0f5970e9354b8f	doi:10.57745/XXTETJ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721309
721234	Serie_10434.txt	13197411	c90c74ceecbdd01acccd8f0e4d89a4a9	doi:10.57745/DNU7P9	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721234
721068	Serie_10435.txt	13209781	bd2b42ff0eb068381008cba8657bb746	doi:10.57745/RD2CEB	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721068
721164	Serie_10436.txt	13177000	1c9f0ce2345f55b522d3c4844370e31e	doi:10.57745/EQGSHB	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721164
721194	Serie_10437.txt	13186255	746cf7ae3c03e9272a91fb4f72600dc8	doi:10.57745/W4PM4K	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721194
721251	Serie_10438.txt	13195473	289123ca6fb1dfb506709e864d110f0e	doi:10.57745/DKD2HK	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721251
720621	Serie_10439.txt	13220848	4b2b0a7d98d5e7ceda6345eeab1d6622	doi:10.57745/VJ9XFN	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720621
720923	Serie_10440.txt	13224552	dc7985d46442fbc80093b0daeb9d04e5	doi:10.57745/MWCOLD	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720923
721178	Serie_10441.txt	13211848	25e106f418fa61806a1aa27f9b3b3690	doi:10.57745/TIAIVY	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721178
720701	Serie_10442.txt	13232983	e2486dc006fc6ecff8afc32182f1bbcf	doi:10.57745/QFZYY2	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720701
720518	Serie_10443.txt	13199682	930d1457d7998bdd9444b45620a79cf5	doi:10.57745/4SZAMU	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720518
721279	Serie_10444.txt	13217933	5779405a15d60f1ec0726215f1224dde	doi:10.57745/II5A7D	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721279
721192	Serie_10445.txt	13209233	535b443973e511217bbdf5d66eb3792a	doi:10.57745/6FY8NI	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721192
721314	Serie_10446.txt	13209573	3d872bed89ee3c76397dc5a43e582a13	doi:10.57745/3XIVYS	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721314
721198	Serie_10447.txt	13227180	39ce7fc736b3d127ea69bcf9255cf417	doi:10.57745/X9HRNQ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721198
721180	Serie_10448.txt	13204863	602a9030e11301a581a339978c3df39f	doi:10.57745/3GVRC9	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721180
721203	Serie_10449.txt	13190921	84f21d6d33ddafaa33f0d1a2bfdc6983	doi:10.57745/4OSLRZ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721203
720804	Serie_10450.txt	13228341	cedfd891c40d947059cbde594a9d03e7	doi:10.57745/MW6SK7	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720804
720824	Serie_10451.txt	13194311	dbe35058d43b4e01d49460a729c81cfa	doi:10.57745/NXXJ7S	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720824
721080	Serie_10452.txt	13186079	2fa9136315e00e4894a398ba403bfbd7	doi:10.57745/OOLXGK	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721080
720947	Serie_10453.txt	13207071	810b3f72cd5808303a557c000210f6b5	doi:10.57745/EU6EK1	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720947
720614	Serie_10454.txt	13164461	485fbab6a8bba8e203e6600293b30484	doi:10.57745/IB7ONI	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720614
720549	Serie_10455.txt	13190483	2467eae186e262c6a94d8530052f4540	doi:10.57745/EBGXKD	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720549
720645	Serie_10456.txt	13188358	1bc29357a298fa2054689fab871d5224	doi:10.57745/HK3KBF	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720645
720858	Serie_10457.txt	13184022	d5a8dde70be6ba10d6cb3b22902eb99a	doi:10.57745/LBZGGT	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720858
720920	Serie_10458.txt	13187436	8c792567b596aa447345007ece24b0b8	doi:10.57745/XCX1TQ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720920
720708	Serie_10459.txt	13210055	19b0cd6c34f842612305a90b8c96fa85	doi:10.57745/YJRKU6	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720708
721318	Serie_10460.txt	13203486	583bdfbd823a4e5f620ad6d4972a10f8	doi:10.57745/E2BXF5	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721318
720410	Serie_10461.txt	13212406	da236421dc45647ff75b54fb8f5532c3	doi:10.57745/5020K7	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720410
720508	Serie_10462.txt	13188751	7985175898e08dcf633a81cedab4521e	doi:10.57745/IOBLCO	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720508
720933	Serie_10463.txt	13209379	7b08f98a98c95199ad41413aa9989339	doi:10.57745/DXV1KU	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720933
721090	Serie_10464.txt	13209685	cf22a7fe4e66084bf245cccfa690ec99	doi:10.57745/T4Y0SU	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721090
721166	Serie_10465.txt	13186336	f427d2b511547089f7487b9efc881357	doi:10.57745/CRW7MA	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721166
720992	Serie_10466.txt	13190233	d0d908b77f2d7080578de8506024f0f2	doi:10.57745/EZRMQQ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720992
720430	Serie_10467.txt	13164981	75085f9dbb6d66a576c9bd33af2aac60	doi:10.57745/JQ6VY1	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720430
721159	Serie_10468.txt	13184832	119ce01f40eb7dcf9e0967fc239df109	doi:10.57745/UDQ61X	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721159
721223	Serie_10469.txt	13185408	09d2140220716bce2eabdc000ff20041	doi:10.57745/1FSZCZ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721223
721123	Serie_10470.txt	13235708	8ac968021fef5f13b1a1cfb1c2a5fa45	doi:10.57745/KYLBFO	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721123
721049	Serie_10471.txt	13225019	7ce7567c25c9fa08e3644cf5202bbc3d	doi:10.57745/5B8ET3	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721049
720652	Serie_10472.txt	13206120	7ff58476c5d9f97a4fc88c564ee6fa05	doi:10.57745/UMOTUL	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720652
720978	Serie_10473.txt	13197115	21b4d6a757c42f4dd7df775837025ad2	doi:10.57745/E69DDK	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720978
720512	Serie_10474.txt	13178370	94baf8925a1a5e7fa5d1494b303b85bd	doi:10.57745/QCGW6M	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720512
720840	Serie_10475.txt	13194870	bc392b642d7872536dd2feb8c557d4e7	doi:10.57745/H50PKB	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720840
721313	Serie_10476.txt	13192243	7da42f9ebdae3761ce0500e810aa41fa	doi:10.57745/YGJEIB	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721313
720969	Serie_10477.txt	13208898	150f5161a68dc61430662542f6320fc7	doi:10.57745/PSOBJV	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720969
721200	Serie_10478.txt	13190468	7e6205fd7707197a0cb88ab003b4e661	doi:10.57745/UGU4SR	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721200
720981	Serie_10479.txt	13228672	d0ed5fad3e9a32949ee82d12b536c142	doi:10.57745/XLGQIW	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720981
721269	Serie_10480.txt	13182110	18182952d4bd9238a4f27e153e965a18	doi:10.57745/EWOQB8	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721269
721311	Serie_10481.txt	13215833	f0eee0525e43b16d1ad048d5b87283e1	doi:10.57745/PK8Q4Q	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721311
720889	Serie_10482.txt	13219416	b8a66832a866005707064868ca827a8c	doi:10.57745/4EFI2E	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720889
720418	Serie_10483.txt	13197547	440f352d72cce80d4fcbc35f93721354	doi:10.57745/JBZBYW	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720418
721139	Serie_10484.txt	13207552	4bdb4037321e43584d3e7c3107b0323b	doi:10.57745/K8LJOR	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721139
720641	Serie_10485.txt	13215225	f76e4dd4b726309045b31fd231da5b0d	doi:10.57745/JPAQHJ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720641
720594	Serie_10486.txt	13218229	4891e4e7d0589f4f5c06e167ada4197a	doi:10.57745/DNYF69	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720594
720925	Serie_10487.txt	13224871	cef6e675793885cdc6f9e1b0443f4ea1	doi:10.57745/IOQFWG	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720925
720859	Serie_10488.txt	13165345	6d162b074d40f10cc6ecb336f785563c	doi:10.57745/HGDB4Y	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720859
721042	Serie_10489.txt	13191419	27909988ab07642215582fccbf6d49e9	doi:10.57745/PN692W	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721042
720568	Serie_10490.txt	13209937	ca74b162328dd06e6b698f3b8d8af2d5	doi:10.57745/0EM2U0	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720568
721201	Serie_10491.txt	13177483	73de6733c7ec07ed74b4cf7bda6496ee	doi:10.57745/9PVMQQ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721201
721084	Serie_10492.txt	13230381	43d5544fd71c2d26725941d56cd62a8c	doi:10.57745/UB53CD	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721084
720787	Serie_10493.txt	13203131	580392eb23ccb2d695ccf960358e3fac	doi:10.57745/AR2UC0	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720787
720562	Serie_10494.txt	13186765	b98336f97eb5791e5dc82a0bb4ae6ea0	doi:10.57745/1NIVOO	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720562
720545	Serie_10495.txt	13195029	8e73c6e0555f7fd126cfe4ab3101d740	doi:10.57745/LOHN3X	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720545
721031	Serie_10496.txt	13183189	c8c3e194993d98297fbf3c19d270c269	doi:10.57745/NUFWWJ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721031
720521	Serie_10497.txt	13195301	3b0c548fc88292e56b00fc24c75dff3c	doi:10.57745/AUJIJP	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720521
721148	Serie_10498.txt	13212317	b02e04521e9d9d6e08cc5bbc619ed445	doi:10.57745/U9IWJZ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721148
721202	Serie_10499.txt	13216179	6b5697e43ee4c0df2e76a37f4df22105	doi:10.57745/FRFTNF	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721202
720862	Serie_10500.txt	13221466	ca6fba7773f3555aafcdd0757df992f8	doi:10.57745/TE2GZ5	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720862
721293	Serie_10501.txt	13215496	7d2070ef3ae463ac243c6d3201aac91b	doi:10.57745/G3BBZ8	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721293
721169	Serie_10502.txt	13203808	1282a47296735d689c835ae6dd4641f7	doi:10.57745/1CJBLO	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721169
720480	Serie_10503.txt	13197786	68350ac3146750ab2b14093e4c9c3839	doi:10.57745/IYWTMP	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720480
721273	Serie_10504.txt	13188188	2caf823b0407198d3e2c72e56d36626d	doi:10.57745/WKKT4L	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721273
721260	Serie_10505.txt	13243821	1fac201dcdf6ef80303217e667b1da13	doi:10.57745/3G34KM	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721260
721291	Serie_10506.txt	13195237	6e7e3aff89e839c52e07701c3d3a92ba	doi:10.57745/M949ZA	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721291
720525	Serie_10507.txt	13205158	1893acb585fbb253c853e951752c5dec	doi:10.57745/ECGRBJ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720525
720443	Serie_10508.txt	13210867	2a1e25609427f30912111bc1ee5c58bc	doi:10.57745/APGLTF	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720443
720830	Serie_10509.txt	13192599	1fee64340f0e36c49e5e4322aed18a4c	doi:10.57745/CYIQW0	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720830
720786	Serie_10510.txt	13198204	b3d23bf70cc2f1ecd7e0da3b6eea5b7a	doi:10.57745/PYB0K5	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720786
720704	Serie_10511.txt	13211094	3b8c1405a900ff6a0884f27ee19d4f83	doi:10.57745/B9LU28	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720704
720937	Serie_10512.txt	13234641	c63329ea8d03d5312ca29a61f75e65cb	doi:10.57745/PNMOKH	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720937
720897	Serie_10513.txt	13210854	94d8cb4116fc8d5eba0630c39a9608f1	doi:10.57745/M9BPRK	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720897
720503	Serie_10514.txt	13203821	e4631f7676f5e90a26a440fdf6842741	doi:10.57745/OUQJSE	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720503
721130	Serie_10515.txt	13199037	da63d1d35645acb31188b65e440b28c2	doi:10.57745/GQPEKK	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721130
720726	Serie_10516.txt	13194884	39126386d68844bb5164bcb96abe4de3	doi:10.57745/MJB7ON	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720726
721146	Serie_10517.txt	13217940	2e38def7e9b81bb55737e44ccd4bb904	doi:10.57745/WU1AUH	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721146
720935	Serie_10518.txt	13205609	a95a762304b282437f266121fdf99cd1	doi:10.57745/PR7IRL	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720935
720468	Serie_10519.txt	13202993	02bfcb07ae9a6f94bd90b3dc5f52e929	doi:10.57745/5LBCAP	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720468
720422	Serie_10520.txt	13236037	b231dca08f59eebbcce25d0fddeb7b2a	doi:10.57745/VXERTW	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720422
720649	Serie_10521.txt	13221511	22a1867d0fa78ff1833a4f10baaa2d98	doi:10.57745/C83ARN	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720649
721190	Serie_10522.txt	13191015	a59a155b33c4596a3afca718238072ce	doi:10.57745/SPICEY	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721190
720942	Serie_10523.txt	13202623	a573dd5651be471d47593622fddef4fb	doi:10.57745/02ACYS	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720942
721226	Serie_10524.txt	13187988	69cd3516b8573230eea747c25dad7771	doi:10.57745/9NPSUM	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721226
720510	Serie_10525.txt	13214163	b79b446788e4c68764b37d4e97458168	doi:10.57745/VLDJSB	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720510
720425	Serie_10526.txt	13234166	a51317da7297f060fb24d15e7749bf34	doi:10.57745/7LAY7J	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720425
720915	Serie_10527.txt	13203844	e887c61cf6429e836f75532ddcd65832	doi:10.57745/YUZBHT	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720915
720833	Serie_10528.txt	13194075	fbea3b56d1e94ae3f624b32c9bbe72fd	doi:10.57745/1XZJX8	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720833
721183	Serie_10529.txt	13220172	ac263701ce78f429bf8a5ed22cc53aea	doi:10.57745/VRZPYF	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721183
721085	Serie_10530.txt	13176121	ab22b9079a9d0eb9c89e7245254e93db	doi:10.57745/VURZRJ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721085
720949	Serie_10531.txt	13197553	5c96481ff7f1cde3730efe10432725a2	doi:10.57745/ZQNGJ4	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720949
720985	Serie_10532.txt	13199100	5402cf6095333f1744d62aea9f45ba23	doi:10.57745/ME4SAZ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720985
721110	Serie_10533.txt	13195967	059a64b2ef8f5671fe6388e2d1e59b18	doi:10.57745/MLCV3R	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721110
720702	Serie_10534.txt	13175329	94e889f1d02f44d5083a986d726d0b9b	doi:10.57745/VBQWUY	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720702
721245	Serie_10535.txt	13206711	bae4adea70380e281969666bf35b2856	doi:10.57745/CHJQBY	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721245
720861	Serie_10536.txt	13208100	f2c847674f76f39bcfb2bbaa2da5c744	doi:10.57745/S28BDS	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720861
721356	Serie_10537.txt	13197135	0b8c7f127f1e1bce125fe2172b0be268	doi:10.57745/YLFBJT	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721356
720873	Serie_10538.txt	13213504	8dc9900b93fd17e6278418125160fb95	doi:10.57745/TKTJC0	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720873
721119	Serie_10539.txt	13227011	c9e94c0d1a0f0376c065e2573624a12e	doi:10.57745/EGOTQW	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721119
720856	Serie_10540.txt	13202115	db0b1dd28f85995ba23434581c3b159f	doi:10.57745/26IQVH	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720856
720959	Serie_10541.txt	13198465	7f963880c2cc6dfce0b42284f3a9ca75	doi:10.57745/GMINV1	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720959
720402	Serie_10542.txt	13202147	d5ba31aabbea94d90a350804cbf6792c	doi:10.57745/V0GCW2	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720402
720542	Serie_10543.txt	13196902	956d3a981227dd744917c73eb97164b6	doi:10.57745/1IMXZ1	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720542
721227	Serie_10544.txt	13203224	d731e7816084d24f5b9cd90f2495d54a	doi:10.57745/L3SC7F	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721227
721167	Serie_10545.txt	13229986	26ea3beec316bab5f20f1ee186e7c698	doi:10.57745/TEAY0K	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721167
720683	Serie_10546.txt	13201962	915e5d9432938fd7b5207f96051c9c3d	doi:10.57745/MKLBI1	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720683
721372	Serie_10547.txt	13208261	8b38f7daf9ed3a88c69c0bf68a48d606	doi:10.57745/ZM3AKJ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721372
721140	Serie_10548.txt	13187197	bd49850dfa52879b6ebd54ddfd1875a2	doi:10.57745/X3R435	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721140
720615	Serie_10549.txt	13175057	b95a9af7a3e47750c5e774a9c3841e43	doi:10.57745/TCFAWG	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720615
720428	Serie_10550.txt	13212765	46803391ba3a14564405c52e6e019daf	doi:10.57745/7OKNRB	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720428
721288	Serie_10551.txt	13217811	13e8b0789a1eb35fef42e2067a3fe4d7	doi:10.57745/MR8NBU	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721288
720499	Serie_10552.txt	13211967	916a6bfb2cc692c87e55232cc9f706b4	doi:10.57745/PCC0XG	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720499
720965	Serie_10553.txt	13224335	c3136226103eb398e3dad1dee6495b62	doi:10.57745/ZXAXNO	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720965
720541	Serie_10554.txt	13199517	7562b22a6c009ca42d132823077e2e42	doi:10.57745/OTKRRP	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720541
720958	Serie_10555.txt	13196007	6f9c3abf442c0127e10a7ce6e31acef6	doi:10.57745/93XOIL	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720958
720763	Serie_10556.txt	13183953	acee306afd7d8fca6fe12e8c33bcb2d2	doi:10.57745/GAFYWW	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720763
720563	Serie_10557.txt	13211117	eebda3a0f3b14c5e64380dd98226a562	doi:10.57745/GZWOHH	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720563
721051	Serie_10558.txt	13164139	bdcc9c9161ac9975a110e2e9e62413ea	doi:10.57745/UYGNN8	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721051
720736	Serie_10559.txt	13164503	1891e54e204b84c4ea39e82f14571aed	doi:10.57745/PAYUNO	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720736
721217	Serie_10560.txt	13207052	066c78da4518ff56acc871dee8b65d62	doi:10.57745/KIHYJJ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721217
721018	Serie_10561.txt	13194074	4f93cb83fc0777e148ffab5c0a65a407	doi:10.57745/LB3Q0R	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721018
721045	Serie_10562.txt	13226163	9236bd9a677a25472eb7a0afcafe1a81	doi:10.57745/UQPYMK	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721045
721256	Serie_10563.txt	13230514	7b70f81eba933eddbf837547537f0a45	doi:10.57745/RC3CVU	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721256
721172	Serie_10564.txt	13211192	f3d101f173a1cd478aa949f960cd3a34	doi:10.57745/QBDCGW	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721172
721384	Serie_10565.txt	13211660	faa02b128b0ad421c61912541d1b9f72	doi:10.57745/9ADNYZ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721384
720574	Serie_10566.txt	13169466	8f172da013cf06be6ea93474647d29a0	doi:10.57745/L50V8W	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720574
720619	Serie_10567.txt	13212316	2da5b8b2143b43853f0ee813da788fc0	doi:10.57745/INZTVO	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720619
720882	Serie_10568.txt	13220052	10b2f2920610ca4d1774c9d24f866569	doi:10.57745/5ADLCW	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720882
721282	Serie_10569.txt	13230164	6941affac6f3c1f7a120d981ac8a3108	doi:10.57745/GYCVDW	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721282
720635	Serie_10570.txt	13203319	8c7ec4b9bb911d0abfd31486157416ad	doi:10.57745/6D6PFQ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720635
721325	Serie_10571.txt	13228094	5c0c76520ac483b130cf03dd40d54ede	doi:10.57745/FRVMLL	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721325
720590	Serie_10572.txt	13221447	4ed66caf890858792c2aef7afa25386d	doi:10.57745/AV7PVP	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720590
720752	Serie_10573.txt	13189789	b52c17a0b95c3f967b7dca838fd35bfe	doi:10.57745/QLANEY	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720752
721229	Serie_10574.txt	13179944	296d680d4563335a86aa79a875b6ef0d	doi:10.57745/A7XXSY	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721229
720688	Serie_10575.txt	13247783	8f35923a93e30be6c08465547cedbabc	doi:10.57745/UN7ELU	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720688
720538	Serie_10576.txt	13212231	f0664340da60b97173cb293d014c5142	doi:10.57745/YXFO3J	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720538
720442	Serie_10577.txt	13214056	d2d931b2424c088999ed86058a0142ba	doi:10.57745/KNWKYJ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720442
720823	Serie_10578.txt	13215880	087033a3390ae9a322b805c7a8a5158d	doi:10.57745/DQCOFH	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720823
721350	Serie_10579.txt	13199234	dbc0a76206cfb9a3b4c36f676cfd21f8	doi:10.57745/GRB5MZ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721350
721228	Serie_10580.txt	13212166	969f8ab365dddbe8c4d7feb9ab78a260	doi:10.57745/UGI04E	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721228
720515	Serie_10581.txt	13206952	d76c380edc79e0c5b3325c8462d0e3df	doi:10.57745/LW4WMZ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720515
720691	Serie_10582.txt	13195400	1c13c11a863236ab0cdae903130fb6b9	doi:10.57745/TTP3GG	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720691
720551	Serie_10583.txt	13229804	3a5ce7287e378b6980bfe632d12df13d	doi:10.57745/BDCUXV	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720551
720678	Serie_10584.txt	13206627	686a72689b0109bbd2b363662bf6a797	doi:10.57745/S4AK9B	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720678
720789	Serie_10585.txt	13210579	78d4f929f081c152f4f28ba5d53eb795	doi:10.57745/D9V3GS	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720789
721193	Serie_10586.txt	13188003	1a878e384977e125e4cebc067ecf97c5	doi:10.57745/XTNBGB	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721193
721361	Serie_10587.txt	13192501	63c0d3a7df9299a7ab0e33ffbddcfade	doi:10.57745/5LEMH1	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721361
720956	Serie_10588.txt	13168567	c7b9fa0fed24d77bb3b3cbadb007d816	doi:10.57745/NHN7P9	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720956
720805	Serie_10589.txt	13177280	38936a1f0ae3fa827a9e06e2f344da25	doi:10.57745/S4OX81	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720805
721342	Serie_10590.txt	13224887	3b56b8d437c61e26e4e4a4aee6e76ee9	doi:10.57745/1OHX61	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721342
720661	Serie_10591.txt	13202087	301825a6e1bde28c3cea3779d2c94a01	doi:10.57745/BCEXKI	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720661
721280	Serie_10592.txt	13202063	23f9f825a93372b8c146a88a7f9f5d6e	doi:10.57745/YNJZKR	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721280
721067	Serie_10593.txt	13218971	503750cbcb6ddb93c551346680660a46	doi:10.57745/BM4JAX	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721067
721205	Serie_10594.txt	13196587	42b3f532341b8e21e75d2cff3cb4487e	doi:10.57745/YVVVGQ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721205
720721	Serie_10595.txt	13215379	8b038593fe5ddf06f3a31ed2a88c8f34	doi:10.57745/JVSDQ7	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720721
720714	Serie_10596.txt	13223350	231d2dabff47af53c9b9206984786d2d	doi:10.57745/3ILOPB	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720714
721321	Serie_10597.txt	13177501	6d428eae5421c7e26d33ce9a1b47f8ea	doi:10.57745/RXWY2B	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721321
720791	Serie_10598.txt	13194220	ecda5c2cb72704745101bdb8d0bd1176	doi:10.57745/VWNQO8	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720791
721210	Serie_10599.txt	13213520	0345e6e6850d701ff89a9c00e8a88703	doi:10.57745/RP2JIX	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721210
721145	Serie_10600.txt	13218950	465fc60b813571c6e56e1a5dfdc4f75f	doi:10.57745/QY8UFB	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721145
720876	Serie_10601.txt	13229717	c4b21c1b4c76beccd6a9c0a2b9a480a4	doi:10.57745/AR6CA4	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720876
720818	Serie_10602.txt	13196172	dfae5189f04d82527d30b56ca3567b0a	doi:10.57745/LEYR7G	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720818
721057	Serie_10603.txt	13219782	b6da7a8520231bb88e1e73ba175f8a81	doi:10.57745/OIHTCZ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721057
720829	Serie_10604.txt	13199415	47347fa45883554874398d97fe55f905	doi:10.57745/EDE0OY	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720829
720437	Serie_10605.txt	13187800	d968d0f990c68181310ade03b0a215fb	doi:10.57745/UUKET6	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720437
720735	Serie_10606.txt	13181727	cd457496fac7dcb9b82ef151796b97ff	doi:10.57745/P3E261	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720735
721199	Serie_10607.txt	13202366	5b8f6128a40dbf2dacac118661fac0d8	doi:10.57745/A1UJRE	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721199
721255	Serie_10608.txt	13217036	8dd196651fd6af8f20bc022dc40b20ec	doi:10.57745/KUSNI7	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721255
720976	Serie_10609.txt	13210492	62424c4589cde4a96d2cea1f6b2a51f8	doi:10.57745/HKOPII	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720976
720815	Serie_10610.txt	13195654	b798fec756322b6473630f7f7d9d431a	doi:10.57745/RE0T4Z	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720815
720413	Serie_10611.txt	13196116	be59a57ef4759a120522c060e3436a15	doi:10.57745/7G1ABB	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720413
720618	Serie_10612.txt	13174326	f3ad94f8aa8c52d75619e91af176abf5	doi:10.57745/KH99RE	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720618
720531	Serie_10613.txt	13210369	a3b2e5ffa8555aa298c9cb35b04f472f	doi:10.57745/RKXU4C	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720531
721156	Serie_10614.txt	13204912	5b5d6e69177184a05881b7881674d789	doi:10.57745/JNHOIW	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721156
720500	Serie_10615.txt	13190408	7358d37cc98787894a3baa40d79abe46	doi:10.57745/BDTBNO	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720500
720580	Serie_10616.txt	13207995	47ff8468f87bc69eca82cf5147c6f381	doi:10.57745/B8L6NG	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720580
721351	Serie_10617.txt	13176035	30bc0a645c2d4dafb8bcf5a43fbef23b	doi:10.57745/8CWC21	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721351
720828	Serie_10618.txt	13176527	cd1a0ff1e82bd3197a7ab3bf2f46372e	doi:10.57745/ICYDRQ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720828
720699	Serie_10619.txt	13202445	6cc118abdab761b73697a03f071dea62	doi:10.57745/FEIWSZ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720699
720809	Serie_10620.txt	13182547	254dde966c085afe98cfbb05c274f0f8	doi:10.57745/KXNOWO	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720809
721150	Serie_10621.txt	13197351	861da47975273b68a7e63364747cba7d	doi:10.57745/QBEBTY	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721150
720757	Serie_10622.txt	13186898	25b27f00617ea4c85bb48a85779a06c9	doi:10.57745/QAROAN	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720757
720761	Serie_10623.txt	13208610	4bde2c83a13f79ef12ee43378ecf5495	doi:10.57745/N3SWLG	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720761
720916	Serie_10624.txt	13216019	7a9adc2ad44ac21ed2f8b79e59a831e9	doi:10.57745/DZQGQJ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720916
720927	Serie_10625.txt	13207533	a00b79ed662fbd2f9959dbf9b074aeaf	doi:10.57745/KXQEIR	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720927
720997	Serie_10626.txt	13230077	29edcfa71c178e17fcc6962af1adbcd7	doi:10.57745/BZFMGA	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720997
720822	Serie_10627.txt	13209646	7d58828f52a5886a1c51ebf8d5b375e8	doi:10.57745/ZPDM7Q	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720822
721349	Serie_10628.txt	13195064	8a6f77d338c88f5bde82816db3f1e8b4	doi:10.57745/SK39MN	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721349
721022	Serie_10629.txt	13208556	0c15a9f534ad02ba61bd4c3abe896aec	doi:10.57745/SXATEH	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721022
720554	Serie_10630.txt	13206658	e0d4b594b90345d45be525bf96105c9b	doi:10.57745/LJJTDM	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720554
721286	Serie_10631.txt	13222486	b1a503d4006f909d5ed167ca7de3d361	doi:10.57745/QNB5I8	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721286
721315	Serie_10632.txt	13219881	c3c991b6c6d94c4e9821ab0154a78310	doi:10.57745/JII2YK	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721315
720526	Serie_10633.txt	13227510	594e97d2c62d26fed24b844b15aefb29	doi:10.57745/XNWYAW	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720526
721160	Serie_10634.txt	13229657	b4188aa57f68cbf8501e1dd4f1e51142	doi:10.57745/GKJOPI	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721160
720463	Serie_10635.txt	13198073	5b688549131b5efb7e4554d9df7d6c08	doi:10.57745/FJVKWZ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720463
721242	Serie_10636.txt	13195955	49448e134416a44460a90197f34ccfc9	doi:10.57745/MVY8PF	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721242
721007	Serie_10637.txt	13215602	9f9b6c1438e3afd1c167cc0eb211e834	doi:10.57745/IXWSOO	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721007
721188	Serie_10638.txt	13193332	b7c36e4097a38afb90a5ba9aabe5fa3b	doi:10.57745/HMLG0F	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721188
720854	Serie_10639.txt	13196501	6c6b794bb3071f1051d0c4b48e3730aa	doi:10.57745/9BGRTV	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720854
720707	Serie_10640.txt	13239371	e9561a8c9c9c27c58649945dfcc096b1	doi:10.57745/3AXE9H	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720707
721107	Serie_10641.txt	13222759	6acacbfb14c16c3edf8055115dcf71bb	doi:10.57745/4WNRGI	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721107
720769	Serie_10642.txt	13212458	db4ac7ccfce9b5ab596c663dca17c27d	doi:10.57745/SNXVZQ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720769
720470	Serie_10643.txt	13173610	3a21c06e70c868c76e2f1453030605b4	doi:10.57745/YS5PKP	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720470
720537	Serie_10644.txt	13202666	baed7a3d961a6d8735f5b1be49c39d38	doi:10.57745/JIQOKM	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720537
720474	Serie_10645.txt	13186847	3b26bfed75e85d8146e855a4c560b7c4	doi:10.57745/XAVVTG	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720474
721388	Serie_10646.txt	13224027	e2923815b460e484b74d976347b88339	doi:10.57745/KXHZPP	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721388
720967	Serie_10647.txt	13220950	c6c37c41f2b9278d0afab727d3728249	doi:10.57745/JYH6EB	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720967
721151	Serie_10648.txt	13194302	cee304cd7c604e8c5e2ff48220a85141	doi:10.57745/CV8HF8	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721151
720489	Serie_10649.txt	13209990	bff21968f1aefc5e480d1f0036474e23	doi:10.57745/8WSMKA	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720489
721087	Serie_10650.txt	13235805	44c798ee7dd80baa9f1fdae26fe5e4fc	doi:10.57745/EQURJ6	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721087
720571	Serie_10651.txt	13185806	4a8bcb626bb2ea6a978a288fcea51e51	doi:10.57745/TZ8ZP1	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720571
721038	Serie_10652.txt	13214142	103e84e4e51203fb79b70f925302e458	doi:10.57745/E5K4Q5	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721038
720473	Serie_10653.txt	13211885	48a00993d4516d0027fb851097360468	doi:10.57745/2VGCZ7	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720473
721121	Serie_10654.txt	13223536	4de7ccaf731820125130eb5273c09989	doi:10.57745/5AANTH	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721121
720497	Serie_10655.txt	13200462	e9cbd8042c45e9e76544c56554f29fcd	doi:10.57745/ZO98VA	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720497
720419	Serie_10656.txt	13205927	3239a9c2cd5719c9a1cff8c8f92f677f	doi:10.57745/ZH3C8J	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720419
721271	Serie_10657.txt	13191909	25176875e9bd3a0655cb145c3a3b61ec	doi:10.57745/QKSRU7	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721271
720974	Serie_10658.txt	13207314	935261ac0b7bb4d28781206f334a0917	doi:10.57745/FSKBN8	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720974
720555	Serie_10659.txt	13190955	dba097106af331c597d6c89e43ce8636	doi:10.57745/EGB6DR	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720555
721061	Serie_10660.txt	13236799	c251ea6ecc7494c0c80d6dc9b9246e89	doi:10.57745/6KIU0W	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721061
720547	Serie_10661.txt	13214997	fa25de4c5c7f232f41a465d9c16fdc91	doi:10.57745/T8F7CZ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720547
720698	Serie_10662.txt	13190984	2158c8348cfec5a7b7c3f7846029f138	doi:10.57745/CDNRMH	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720698
720513	Serie_10663.txt	13186706	13cdbc361c0a9ff52d4e3e69eee8b74d	doi:10.57745/XXWC3N	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720513
720578	Serie_10664.txt	13220670	8591e93c022ab0d7738fb51d1b6b6d96	doi:10.57745/K5VVOM	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720578
721312	Serie_10665.txt	13186981	610cf6ce040d832ee1290cba8b5df9f2	doi:10.57745/T9KNN9	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721312
720806	Serie_10666.txt	13194615	410e8bb92a0f1bd5b1221c18e79da2f1	doi:10.57745/MR6NMT	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720806
720821	Serie_10667.txt	13168938	f7496732138a5a028daa0c15878c5c8e	doi:10.57745/RLKCBH	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720821
720866	Serie_10668.txt	13160357	ef560998840b81579bf36047a75a65ea	doi:10.57745/8KVZPZ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720866
721137	Serie_10669.txt	13213514	ec60d314025bddd339aad3b82f3c0a2f	doi:10.57745/STAKBF	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721137
721374	Serie_10670.txt	13215766	a9d09bc0e5d8b8da7ae7e14b1fcd25a6	doi:10.57745/3UX1PL	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721374
720745	Serie_10671.txt	13219368	088cd0e61fbb5e2fe3108ea907b07d16	doi:10.57745/1VOM2Z	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720745
720624	Serie_10672.txt	13237149	bb75e8005650d39e30805c85e8a8f796	doi:10.57745/A4FC1K	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720624
720759	Serie_10673.txt	13201020	83f2a723c95a3e9d1170e935bd8c0686	doi:10.57745/TODJTF	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720759
721109	Serie_10674.txt	13219609	8e33b3beba1e16dcd4fca33c3e28fe82	doi:10.57745/FRFDG9	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721109
721155	Serie_10675.txt	13201614	605bd0c53fa0c30732c3b59d289f839b	doi:10.57745/XGIJ2P	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721155
720834	Serie_10676.txt	13207372	e9a944fed1f4dda7a2bfe1d0f2c05043	doi:10.57745/ND5Q4S	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720834
720742	Serie_10677.txt	13188306	67ef95d24e595448bdb6b5c015e502a8	doi:10.57745/SAWK2U	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720742
721308	Serie_10678.txt	13195843	91f7173515a8df99d3951048c32af843	doi:10.57745/5OJSYY	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721308
720712	Serie_10679.txt	13216903	a6fbff972d85eed854aae9f4c0b72750	doi:10.57745/556X85	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720712
721065	Serie_10680.txt	13206778	7c24251f3b6a322c3996bf355009d90c	doi:10.57745/8ZDKAR	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721065
720604	Serie_10681.txt	13222721	bf5764fc08a13eb688a281075f2042e8	doi:10.57745/FCPJGD	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720604
720774	Serie_10682.txt	13195746	a3ca5b17aff8b549945e649831f0d4f4	doi:10.57745/D4RCWD	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720774
721071	Serie_10683.txt	13199349	b28f96b206a3705b336c8253dd3438fb	doi:10.57745/G517IJ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721071
720629	Serie_10684.txt	13192273	c394afb184552ad87f932688e82bb895	doi:10.57745/Q0IZIQ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720629
720754	Serie_10685.txt	13198766	c10f4ad0a618ec787fa40155cf9daf86	doi:10.57745/J8JJDK	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720754
720434	Serie_10686.txt	13209578	7b92d5b258b604756a34bfe93c79e6cc	doi:10.57745/42CLOY	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720434
721174	Serie_10687.txt	13209144	07a25408151010011b54c77b83138234	doi:10.57745/CZ13Y1	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721174
721319	Serie_10688.txt	13209976	8105a899689c366c1516e2e172ac92eb	doi:10.57745/MVQA1N	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721319
721249	Serie_10689.txt	13223143	56a01605477a86f6e5d30f3f89393522	doi:10.57745/OXYZXR	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721249
720603	Serie_10690.txt	13171516	1bb766a3170625cec0311e8dcf964747	doi:10.57745/CIT4TM	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720603
720718	Serie_10691.txt	13208234	996db9d1ad11c76c21e50612bbe7f4f5	doi:10.57745/XHYMZM	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720718
721264	Serie_10692.txt	13198970	afd83d682a3dfe90951afbe227ecaa50	doi:10.57745/OF6MVA	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721264
721320	Serie_10693.txt	13212629	6da791eb8f8c3dbcbf5e1ef818a35743	doi:10.57745/SRKBWX	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721320
721176	Serie_10694.txt	13230805	8b3e0f74310ed1da514a2e0b8de66162	doi:10.57745/LR8LBT	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721176
720844	Serie_10695.txt	13221109	a2a4328ba7fa3f8a34f0e3ca20eefd15	doi:10.57745/DDXZO7	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720844
720814	Serie_10696.txt	13235428	4a418f32af8110cfe1a99547a3de123e	doi:10.57745/0CGVJZ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720814
721058	Serie_10697.txt	13227868	fe5f134901dd1b8dcb6cbb8b60b10714	doi:10.57745/EL32TN	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721058
721195	Serie_10698.txt	13198807	1af66a018e772ad877a2b57d243b9f82	doi:10.57745/LDM5SK	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721195
720725	Serie_10699.txt	13186969	ee6607fd86c844f4395ff888495f869a	doi:10.57745/ZHMWSI	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720725
721389	Serie_10700.txt	13163620	f5c8dab48fa25c763db2699670ca8ac4	doi:10.57745/PEZ8TL	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721389
720898	Serie_10701.txt	13187760	e2e905392513156a051a29381645c606	doi:10.57745/IEMWPO	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720898
720609	Serie_10702.txt	13202726	e986c361d6024687b8dfa12d254ad6d3	doi:10.57745/EDEJSK	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720609
721129	Serie_10703.txt	13188352	068561a1d90987ae7b7e9598914c96fa	doi:10.57745/ZEEUFG	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721129
720744	Serie_10704.txt	13196696	2f96204e43040b39a46d40ada2c4e7bf	doi:10.57745/R4BTU6	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720744
720440	Serie_10705.txt	13201428	b83709824b9830b0d769f4bbb2b2231a	doi:10.57745/EDP5UE	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720440
721332	Serie_10706.txt	13192632	27718a6d15eca35eb05f1419c8ae97a2	doi:10.57745/6TRBNF	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721332
720857	Serie_10707.txt	13187978	0e4266c75ad5dce5dd8b63b27fe07412	doi:10.57745/G2NMZE	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720857
721025	Serie_10708.txt	13218624	2e1969f42695b6c7f95da46029029762	doi:10.57745/UHWIAT	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721025
721266	Serie_10709.txt	13206623	74d975d0306e76db069a8b18e55ac8a6	doi:10.57745/XIA4T6	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721266
721132	Serie_10710.txt	13191919	a48ab3440c75eea7f04802572eb37413	doi:10.57745/JJM9TH	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721132
721246	Serie_10711.txt	13200921	f647cd723e2cbdff69fa75d88872d885	doi:10.57745/TGUBG8	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721246
720663	Serie_10712.txt	13207250	8f0c61ae1b66da20b2876c873f348c49	doi:10.57745/SUQEIM	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720663
720983	Serie_10713.txt	13215798	fb310c392b48b226d038d6b33c531651	doi:10.57745/7JR2FR	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720983
721149	Serie_10714.txt	13218557	0108939c126a519693a19326e1c7e71d	doi:10.57745/WCY7CX	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721149
721079	Serie_10715.txt	13191797	e2c2b58d37b85dede9e37a06b0467d59	doi:10.57745/MUWPR1	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721079
720878	Serie_10716.txt	13175775	e8e89dbb92384122e6ef5aa09354a64a	doi:10.57745/HTNYBL	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720878
720456	Serie_10717.txt	13180338	d01be8e3de5ce95c3018f350f85ab2af	doi:10.57745/0YQFQW	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720456
720492	Serie_10718.txt	13192871	56678c4bd917963005bb6b95bb8176e0	doi:10.57745/U3HVX4	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720492
720679	Serie_10719.txt	13213337	6ca93bce6ab0f4c04fcb43a33e72a1e5	doi:10.57745/EOLRIT	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720679
720998	Serie_10720.txt	13225158	91d9148b936e16cea14bcf7959789fd3	doi:10.57745/8FXQ1S	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720998
720642	Serie_10721.txt	13219314	9c4d7892c27793a35cac887702f4e866	doi:10.57745/43XGA9	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720642
720444	Serie_10722.txt	13190966	83ae6abd5237d9123572e5961c0c8d65	doi:10.57745/NAHQ88	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720444
720794	Serie_10723.txt	13185988	db1c3e3a708ef1044bf5b807b6f9ca8e	doi:10.57745/RZA2XA	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720794
721027	Serie_10724.txt	13188224	9e491416adf0799121583fea6664e523	doi:10.57745/HX8SVV	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721027
721270	Serie_10725.txt	13200796	2f18e871c6b70f606846b82de8b3d094	doi:10.57745/UYIBTJ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721270
720826	Serie_10726.txt	13210824	a5323f655b692882b4af8272e77d4e2d	doi:10.57745/GVTADG	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720826
720611	Serie_10727.txt	13216555	43c854cd74ee04a77f1d6b2fde082f83	doi:10.57745/C60TPR	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720611
720488	Serie_10728.txt	13219392	e1b70f6a20978186761237c39b323312	doi:10.57745/QZ8A52	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720488
720980	Serie_10729.txt	13205045	405961fd21f28e6631b6e45d53dea073	doi:10.57745/XWHOGC	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720980
721144	Serie_10730.txt	13192579	9052b2a0bd4b239ecc11b7c6950978dc	doi:10.57745/B4FTLT	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721144
720485	Serie_10731.txt	13199478	eddfc417d1380c13eb9a34929d9cdfa0	doi:10.57745/0FFMI3	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720485
720522	Serie_10732.txt	13226528	ac803fb63b133f3b7b11357ec778331f	doi:10.57745/13XCOI	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720522
720501	Serie_10733.txt	13161564	1bd666f9dd9914146d542b101f200919	doi:10.57745/LOJOC4	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720501
721317	Serie_10734.txt	13198401	0ab95f49a814b79930b38f13dff487cf	doi:10.57745/T8PLN1	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721317
720768	Serie_10735.txt	13206259	e080988a4a7fed5e42c289d23edf3d77	doi:10.57745/MZKFCZ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720768
720441	Serie_10736.txt	13196197	6464e88399e69a7d0e48fa971ae0f3db	doi:10.57745/TZVAIR	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720441
720600	Serie_10737.txt	13211691	0b38da1e2c55cadc6e4526c0f4eceed0	doi:10.57745/ZH363Z	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720600
720465	Serie_10738.txt	13211344	e592f741e2f345880c59a4733b0256cd	doi:10.57745/RG1FV8	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720465
720606	Serie_10739.txt	13202688	0f41743fac438a277840caee60d0a317	doi:10.57745/AINGT6	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720606
720408	Serie_10740.txt	13190808	95f97976bc2e50f0bf0d33634f42d714	doi:10.57745/ERBYCY	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720408
720491	Serie_10741.txt	13189720	bfb6478cec93f9ed18bcd18a869be2d1	doi:10.57745/QT67ST	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720491
720820	Serie_10742.txt	13192710	478a2908bd807cd22b711034fc84a676	doi:10.57745/TY3EKK	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720820
720680	Serie_10743.txt	13217277	dd8ba57148f85d2540e81fd9a80c5130	doi:10.57745/EB4WIR	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720680
720767	Serie_10744.txt	13196705	02ecdfcfc483b792f7af546055c2cb28	doi:10.57745/JYELEA	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720767
720591	Serie_10745.txt	13198278	d52d42ce49d28f753fcb54b1ea5f045e	doi:10.57745/7QDAUB	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720591
720514	Serie_10746.txt	13186943	fbf434ae0a120b333d49454635ffd202	doi:10.57745/LP1KY6	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720514
720723	Serie_10747.txt	13196898	78658825ac39f05ce9989e59145ecb58	doi:10.57745/KLIBXO	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720723
720749	Serie_10748.txt	13186743	53b41423362dcdb2cb56b7f94de0309c	doi:10.57745/PEBHVI	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720749
720877	Serie_10749.txt	13206662	11f6aa8d43a0acac8da41f6653068d71	doi:10.57745/V6ZAMU	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720877
720546	Serie_10750.txt	13207571	366c2a6202bc7d04ef213328885b8e20	doi:10.57745/XNBXUX	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720546
721093	Serie_10751.txt	13159166	b3b2a1d4100405469d22dd97fc5375db	doi:10.57745/5O16Y3	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721093
720631	Serie_10752.txt	13198696	8369f2e3c6536859327164b0bea32260	doi:10.57745/6FIBKB	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720631
721337	Serie_10753.txt	13204847	c4adad17a6c0d10f3468510ec23f5ab9	doi:10.57745/I1W10I	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721337
720502	Serie_10754.txt	13213459	39b567c3c3218c52a1269d037d7f7823	doi:10.57745/MXZ4B3	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720502
720852	Serie_10755.txt	13212483	3bdc13310ea2ae6b5f333f969d6b6ccd	doi:10.57745/SUKIOB	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720852
720884	Serie_10756.txt	13233650	0ddf8d77258f364d05be9394cd9808f6	doi:10.57745/W4GBIR	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720884
721378	Serie_10757.txt	13194685	eff6fbefbaf808cefb6d3c819a3139b1	doi:10.57745/0ILIRJ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721378
720987	Serie_10758.txt	13171654	07e2ce3ef5870b74a9139b6221391f45	doi:10.57745/IQCIMY	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720987
720561	Serie_10759.txt	13210958	26f3f30f115000564dd15bf282d9618f	doi:10.57745/CMUBHS	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720561
720734	Serie_10760.txt	13198692	1315162b2785ca27b2f9734e6f30305a	doi:10.57745/E67DAF	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720734
720914	Serie_10761.txt	13220296	953a8f839d387e8573cac7bb5624b6dc	doi:10.57745/5D5SLR	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720914
720945	Serie_10762.txt	13224318	5a090c898fb67a51066eb56313e0b8d7	doi:10.57745/9ZQ736	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720945
721219	Serie_10763.txt	13207735	70614709002d9f2ffbc499512d4080eb	doi:10.57745/FJSJDE	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721219
720989	Serie_10764.txt	13221872	64899d9ab30e3ce249343379bacd6166	doi:10.57745/1YLQN1	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720989
720853	Serie_10765.txt	13177602	7f91512da00d8faee414e2ab864185bc	doi:10.57745/A0PV2G	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720853
720796	Serie_10766.txt	13195564	4c2d03037803751632c0168e1999d887	doi:10.57745/1H5HG3	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720796
720938	Serie_10767.txt	13189730	23b9937cd994de6b249546a739f09550	doi:10.57745/TLKCVQ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720938
720932	Serie_10768.txt	13231991	155478d2821a2f206852d78701c5be98	doi:10.57745/44KHHR	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720932
720899	Serie_10769.txt	13193575	7c9f9bd826ebed92cb2c3e3b478d0321	doi:10.57745/XZRYUG	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720899
721366	Serie_10770.txt	13204987	5b885b5c83c214ee33396f6610156ddb	doi:10.57745/VFWAB4	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721366
721002	Serie_10771.txt	13210965	05c72594bfaaf6db5c732a4af69bc867	doi:10.57745/R24UMM	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721002
721024	Serie_10772.txt	13194498	8fc4b54da49c5453eca7fefd27d77964	doi:10.57745/JQY0OV	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721024
720870	Serie_10773.txt	13197132	35ae44b22404ac8aaa04ae6e10dbfc20	doi:10.57745/1EHL4M	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720870
720831	Serie_10774.txt	13187335	1248fe28214865224a97b4e7bcb54073	doi:10.57745/VLFYUS	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720831
721184	Serie_10775.txt	13214322	a0504892ffaafaf80adf40851b92cd91	doi:10.57745/ZNEOGJ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721184
720926	Serie_10776.txt	13214699	402cc39cdb1e9bd6fa19d5dc52531486	doi:10.57745/OAIPY5	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720926
720469	Serie_10777.txt	13192177	7c314028a822a5164de610efa96d8766	doi:10.57745/PBOD4E	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720469
720633	Serie_10778.txt	13205375	8841a003fb5b19738d15dce7c21afc7a	doi:10.57745/ZBE42L	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720633
720539	Serie_10779.txt	13190555	6a61fd79c84a5e17a4762e49f5a622f6	doi:10.57745/F14LQX	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720539
721189	Serie_10780.txt	13194067	54397bba358c8f06083e9c7da8d62be1	doi:10.57745/770QO5	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721189
721360	Serie_10781.txt	13208754	4cf99b1a032c54d348b864b4acbef8b9	doi:10.57745/UAN2DQ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721360
720968	Serie_10782.txt	13175223	5d5383e3e6f74c7a950f4c3bdda69704	doi:10.57745/3REIGH	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720968
720602	Serie_10783.txt	13225980	c58f50f65b563f571e1c9a3ac880c004	doi:10.57745/VUJRPT	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720602
721295	Serie_10784.txt	13193196	70ba3ef3bdf8d68368db0eb6aa0e7cdd	doi:10.57745/YD2D0P	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721295
721287	Serie_10785.txt	13207391	9b612a7a496ec06fe3eddbccc1ecd435	doi:10.57745/VWWUX6	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721287
720684	Serie_10786.txt	13191363	ba840d4e1e1134cf9a34571b2b320b6d	doi:10.57745/SFRMH5	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720684
721208	Serie_10787.txt	13209487	7384dd35056ded0292931c8ddaec6639	doi:10.57745/KBOEE8	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721208
720839	Serie_10788.txt	13197461	f33c0d7c04b3f58ba3d015085c3fb1c6	doi:10.57745/EGNBJ1	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720839
720543	Serie_10789.txt	13205672	54e4831b55455fc4f519ef9e6f4dc59e	doi:10.57745/XQUNEZ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720543
721247	Serie_10790.txt	13208418	98de82cad62f1cee2e51519906cb5289	doi:10.57745/OA0VOK	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721247
720773	Serie_10791.txt	13204363	5f1593049a88767e90e35fe51fba5a60	doi:10.57745/OZAM26	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720773
720948	Serie_10792.txt	13221708	2670e9780b2c11d5e10f3520815616f9	doi:10.57745/RWBSP7	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720948
721021	Serie_10793.txt	13200396	6c56d3d13303eea6f2da77a85c799e4d	doi:10.57745/HDLTX0	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721021
721354	Serie_10794.txt	13202315	6e25e3be0b46473eda013c6c64fbda5b	doi:10.57745/MZS0NG	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721354
721268	Serie_10795.txt	13200606	480719511fa9b1d09c99cc4b0c8b3ca4	doi:10.57745/GRZVR3	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721268
721052	Serie_10796.txt	13222202	9550bbffced2b193d791d094a65af558	doi:10.57745/4PNSXM	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721052
721054	Serie_10797.txt	13218083	53a9045a8c4b9a85374f7dfbf0290478	doi:10.57745/3DP9Y0	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721054
720524	Serie_10798.txt	13200883	d88eb536c1c615a8ca1a6f6e43121537	doi:10.57745/JIMLV9	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720524
720893	Serie_10799.txt	13175233	7621e7a904dcc99ca37c75226c48b764	doi:10.57745/QYG9JH	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720893
720705	Serie_10800.txt	13167101	c83c53b6ee70daa6483e777070d9839e	doi:10.57745/6JA5LK	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720705
720883	Serie_10801.txt	13163649	9b061723e392afea6e2bc2c19eb6f994	doi:10.57745/AHC7VB	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720883
720636	Serie_10802.txt	13221414	50b555eafc865d93d8a805b64ccf99d1	doi:10.57745/IBDWO0	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720636
720569	Serie_10803.txt	13203191	797d9163688adf03532eba05334f94ae	doi:10.57745/VVSNXQ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720569
720666	Serie_10804.txt	13229891	b99d40e7874106e189ba35a0820adc4f	doi:10.57745/NIGG1P	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720666
721304	Serie_10805.txt	13226752	9e559b3fe963bb476ca0f4deeaa5e94c	doi:10.57745/PGPLVJ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721304
720651	Serie_10806.txt	13180376	c99fa6367dfcd501f9350f055402a958	doi:10.57745/TX4SH4	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720651
720748	Serie_10807.txt	13178881	e8041da6ba9312f9a3be9c9de768180c	doi:10.57745/CYQE6M	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720748
720803	Serie_10808.txt	13181411	9361dc7dfff30934f5c5882f65075abb	doi:10.57745/7BH1GD	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720803
720946	Serie_10809.txt	13211446	55e5b917c491d33cd0bac58baa337446	doi:10.57745/FUNIPU	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720946
721117	Serie_10810.txt	13204140	e7735fa843b6b2d1a85bb9dead9d77e2	doi:10.57745/P9B0EB	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721117
721365	Serie_10811.txt	13218422	24c5a608f95d89d2b8d99f22e94c6648	doi:10.57745/EOHLI0	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721365
720955	Serie_10812.txt	13206681	1bf02c82dd434ceb72093e1d7327b57e	doi:10.57745/DO8ICO	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720955
720587	Serie_10813.txt	13223895	0be88841cd0c78ea6b07d803d18c2c40	doi:10.57745/YOZ6MR	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720587
720424	Serie_10814.txt	13237455	b2776c75383f6a5ac6c28c2f611cf487	doi:10.57745/ZYFLV5	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720424
720963	Serie_10815.txt	13211686	0c58d4be989da5eada788c0481376b2d	doi:10.57745/4XHS2T	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720963
720690	Serie_10816.txt	13211095	feb8306883270ac32b70ae175661492f	doi:10.57745/RU0DHC	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720690
720610	Serie_10817.txt	13208860	2db4521cebc6fec47766c9a368fb8432	doi:10.57745/XGPDIT	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720610
720461	Serie_10818.txt	13170270	38a376bb4e574db37c00cf74cefd9691	doi:10.57745/FSKH4H	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720461
721265	Serie_10819.txt	13211115	7b034f64cdd4e4d0cebcc9a5cc671d5a	doi:10.57745/NITFMV	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721265
720894	Serie_10820.txt	13189082	98a490b01b6ab2c82800f636b0c1d34e	doi:10.57745/MMHVHG	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720894
720817	Serie_10821.txt	13196672	096f8f2b755ef5338576c4b61f24bcc0	doi:10.57745/YAC0RP	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720817
720449	Serie_10822.txt	13220999	159a5aa925c41d4dfe4c6249074be68b	doi:10.57745/OOFFBY	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720449
721177	Serie_10823.txt	13203519	d06bee209de4a2d6041dde5fb17b3002	doi:10.57745/KSEM4S	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721177
721300	Serie_10824.txt	13209075	f62125fb684aa4e87573bfd777047db5	doi:10.57745/UHGNNP	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721300
720943	Serie_10825.txt	13198201	7ee908760cb27d36fac7b45f3096f7f3	doi:10.57745/RFRNAI	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720943
721290	Serie_10826.txt	13208657	358eb4612f8ec05b01377ea29aac4ef4	doi:10.57745/7Y2COE	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721290
720593	Serie_10827.txt	13170358	7403ecc0f4cb0de1afb5d4e203a2b475	doi:10.57745/C0RX20	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720593
720812	Serie_10828.txt	13206893	eac95d2a9e8207d1e798fefbc34e3533	doi:10.57745/EBABW4	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720812
720671	Serie_10829.txt	13200712	e6ab6483e614a84dbe8d73b0234605e0	doi:10.57745/FWUSBD	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720671
721050	Serie_10830.txt	13193779	86a9de514a1f2948a1b07654bfa93bfc	doi:10.57745/SWIB05	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721050
720952	Serie_10831.txt	13236261	8b1aa2bd4745fa0ed919a66bc1e63fa7	doi:10.57745/Z414HX	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720952
720486	Serie_10832.txt	13222411	0a70e6cdb4007004409418b8fb8f8faa	doi:10.57745/OIOBIX	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720486
721359	Serie_10833.txt	13195248	2c41291ddb3d87f10ca966833f5c7ad1	doi:10.57745/ZQJYJM	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721359
720887	Serie_10834.txt	13197276	6216bf9ff31e40f0fe27bccf54d8abdb	doi:10.57745/ZON498	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720887
721233	Serie_10835.txt	13224968	3776105b1c2c949dad86ea4eba36a832	doi:10.57745/CXAMVG	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721233
721143	Serie_10836.txt	13195944	fad8c76161be5cf9c529d7d282c93851	doi:10.57745/JZS0IE	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721143
721352	Serie_10837.txt	13209892	a5e5da962d757cab7abdec8e343892da	doi:10.57745/NMZIKL	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721352
720588	Serie_10838.txt	13196245	6dd58d0c1cc17967f7aed77f43a5105f	doi:10.57745/JYXEDS	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720588
720842	Serie_10839.txt	13192404	61c4d47d983571c2290d084c482be756	doi:10.57745/X8K9PJ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720842
721035	Serie_10840.txt	13226358	b15588612bc24ad7793b462527c0e2d7	doi:10.57745/QEMLB1	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721035
720835	Serie_10841.txt	13186919	31bdd8d6970814f69fa43a9a5523c785	doi:10.57745/YOQYUL	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720835
721047	Serie_10842.txt	13194239	737d66cb0602794de6e9b6a1f42d4f68	doi:10.57745/YCZREI	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721047
720797	Serie_10843.txt	13156978	4da746f114ab289e0cec88e38386e7e6	doi:10.57745/MJRMRH	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720797
720623	Serie_10844.txt	13225865	69b99d3dffb28787db69dc24ac273929	doi:10.57745/G3CPIC	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720623
720429	Serie_10845.txt	13201559	2e62ed1c5190ed26a4662bf405d0b66d	doi:10.57745/LLJEQB	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720429
720739	Serie_10846.txt	13209475	e9abd9c1233c1489c5e5f98e81d41d6a	doi:10.57745/20YAYA	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720739
720458	Serie_10847.txt	13217324	eb26ab5db94952d114ef0bf233ddf58c	doi:10.57745/EHRLUD	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720458
720733	Serie_10848.txt	13224417	cb33650ef96b4f68b7756954fd5e1c18	doi:10.57745/4AQVUE	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720733
720843	Serie_10849.txt	13212056	a374e89fc2245dc0bdc6ff955071db67	doi:10.57745/KRCA1U	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720843
720667	Serie_10850.txt	13194858	a98f6fa64408a7596448b61d18757f99	doi:10.57745/YXTVW0	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720667
720984	Serie_10851.txt	13194194	b7282e1b59526558955fc6a49fdbf6a6	doi:10.57745/ZF7BNM	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720984
721128	Serie_10852.txt	13218370	89d465e13acc3cf11d3b406c9d44d230	doi:10.57745/SRFHY8	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721128
721094	Serie_10853.txt	13231991	a75936cd8a3e7aec4ff5afc3308ee28b	doi:10.57745/UTEIF5	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721094
720550	Serie_10854.txt	13182907	9875ef18192a9607040acbeb8842ad91	doi:10.57745/O0PYVU	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720550
720482	Serie_10855.txt	13232228	ed9c55f3185db246156b9ba01df9cd86	doi:10.57745/3HFCNO	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720482
720918	Serie_10856.txt	13184392	15bc642c55d1f78835e1afe00cf02bfb	doi:10.57745/VUZ2Y1	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720918
721296	Serie_10857.txt	13194075	6f97cf3aab52e499cc9aa358b4b6e641	doi:10.57745/QQOY4Z	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721296
720548	Serie_10858.txt	13207189	7e6be0519e553f03c7cd795834da9925	doi:10.57745/M1BNAM	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720548
721179	Serie_10859.txt	13165463	4c1c2ec725c3353ecc2c0316a0372ce2	doi:10.57745/EZDQC7	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721179
721010	Serie_10860.txt	13204695	eb385b363ad0a63029dcac7cd537d305	doi:10.57745/HR33SQ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721010
720910	Serie_10861.txt	13199983	e41e15c8c9e7420225ada7a995a3f3e5	doi:10.57745/NNH2CN	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720910
721333	Serie_10862.txt	13226177	16ff6f5fc4fbe182cc089e1b4b8f1e9e	doi:10.57745/VTLOXN	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721333
721298	Serie_10863.txt	13225240	bb9fafdb6102f790118258bb79945228	doi:10.57745/S3Q7CQ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721298
721171	Serie_10864.txt	13218996	4fe91585c175f03a7514c4ec7c80c8a9	doi:10.57745/D42DK2	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721171
721375	Serie_10865.txt	13199373	8f3566d2feb81e74992bba71034bc2e1	doi:10.57745/7SW910	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721375
721026	Serie_10866.txt	13206339	05431649f47ecab1b82ed9ba6dea4c2b	doi:10.57745/X11NRI	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721026
721218	Serie_10867.txt	13175405	d9aa42e064e1d4610a7c60a10bdaac07	doi:10.57745/FTTBRU	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721218
721098	Serie_10868.txt	13212480	a1970bee6cddb7c8626e04c6f2151365	doi:10.57745/6ZVF9G	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721098
721078	Serie_10869.txt	13183205	9f170d133a1a3186c14ceccfada5522b	doi:10.57745/Y9L02L	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721078
721220	Serie_10870.txt	13223422	ad08e0fee70b749ee656b1c1bcc6ab6d	doi:10.57745/S5EA4C	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721220
721275	Serie_10871.txt	13213415	21bc92d44e8fef3030f87456c242d142	doi:10.57745/HTAHDR	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721275
720972	Serie_10872.txt	13211282	749ff40e3a7e9c996a0a4f88f552b000	doi:10.57745/VTDQLR	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720972
721241	Serie_10873.txt	13208999	19e378274866a45017147ba80b4b3d43	doi:10.57745/MLDNQE	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721241
720681	Serie_10874.txt	13186179	04747eaa078b47097ae2fa1730ae3806	doi:10.57745/TRPGLK	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720681
720476	Serie_10875.txt	13190569	9ca8007e0568c2457e480130216f6e07	doi:10.57745/PT5UQA	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720476
720451	Serie_10876.txt	13162946	2b57c705ddb857cd155c30023954bae1	doi:10.57745/61SCI2	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720451
721062	Serie_10877.txt	13200343	fd77469f4bd9f35dd63e0eedf87ede61	doi:10.57745/07C1B3	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721062
720790	Serie_10878.txt	13210255	acea595df1a8a3754b494eea10b16a36	doi:10.57745/O91UZB	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720790
720716	Serie_10879.txt	13202308	fec5f74c4bd5d5833c858747a939c441	doi:10.57745/SJUCLP	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720716
720819	Serie_10880.txt	13211949	ae45ed82af00c043dbe54f7e02f3cd11	doi:10.57745/TTATU5	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720819
721053	Serie_10881.txt	13240366	bb5b339e740cc6a6fc6cb2255dc9a5a3	doi:10.57745/LV5ZSW	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721053
721297	Serie_10882.txt	13174258	6c8a0327b7334b4aa5aea57bef7079ae	doi:10.57745/0NBBXN	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721297
721059	Serie_10883.txt	13182771	87d275e024ca067ba0c5224e2141c36d	doi:10.57745/MRNPPL	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721059
720566	Serie_10884.txt	13213982	47f8ed0f3eaebd63bb9f70206b5a9ef4	doi:10.57745/ONUGKB	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720566
721134	Serie_10885.txt	13219982	79172287651fa51a5336516507b0e1c9	doi:10.57745/KJZPNE	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721134
721257	Serie_10886.txt	13232514	98dad353d14cbdcd891f565c294f31a4	doi:10.57745/AI2XSO	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721257
720426	Serie_10887.txt	13185801	1fdab34163da68cd3fed81e50abfe94f	doi:10.57745/4OGW3W	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720426
721331	Serie_10888.txt	13208537	548d5dd82dff52567fe026407958f49b	doi:10.57745/GJFBLZ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721331
720601	Serie_10889.txt	13196947	b579106094f67157f5c10d1fe0880113	doi:10.57745/TBUE5E	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720601
720608	Serie_10890.txt	13220593	c2e3687e25aadde08b31a6dc4dba119a	doi:10.57745/GOYNRI	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720608
721118	Serie_10891.txt	13202288	97bac7b2f5ff09cee2624eaa688c77ac	doi:10.57745/OWOPQF	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721118
720961	Serie_10892.txt	13207587	20b701af108ddcf5380dafc174fb3680	doi:10.57745/57AB4U	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720961
720658	Serie_10893.txt	13212180	9520ff80bbed30888e127107961d8a5e	doi:10.57745/G51ZBJ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720658
721056	Serie_10894.txt	13173753	4784126385e254ac9346e36f2319d918	doi:10.57745/QRJN1Z	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721056
720564	Serie_10895.txt	13199172	1c6ede65fbe3a2b2bd1eccf8c182c77b	doi:10.57745/ZXI9LE	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720564
720573	Serie_10896.txt	13200412	7cd93c8fc84b0d5b97cb8927d8113dfb	doi:10.57745/KCRYPZ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720573
721046	Serie_10897.txt	13170248	4af220d9a2d165344087329b2da403dd	doi:10.57745/OWCQNJ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721046
720596	Serie_10898.txt	13218248	6b9db575cb52dffdafd9631fb4f25126	doi:10.57745/TSVUOF	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720596
720557	Serie_10899.txt	13195099	2b7d7f60ee26405f1b9c7fcf83ff6568	doi:10.57745/BGOH7D	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720557
720929	Serie_10900.txt	13174138	4c617190b4305d125fd15db1771f3ce5	doi:10.57745/OCFPIW	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720929
721154	Serie_10901.txt	13187532	5788bb7563b3c8b3fcf6948ad5a11ae0	doi:10.57745/XERZCU	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721154
721370	Serie_10903.txt	13214274	8740531a64b53769e94ffa6c5b9be124	doi:10.57745/E0I16U	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721370
721076	Serie_10904.txt	13165772	0b9d1c8fef2a7d594ff2e1757bbdeaef	doi:10.57745/RVJ44N	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721076
721105	Serie_10905.txt	13204260	c0e28b5e80fa90f7ccbac8bebcbd951e	doi:10.57745/H7ZXLJ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721105
720895	Serie_10906.txt	13172866	a2628758c6f7bf1eee04add28032a9e1	doi:10.57745/PDCYCF	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720895
720420	Serie_10907.txt	13193789	f83d10628695b07b8a0d775c7f9bc4e4	doi:10.57745/KQN3QA	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720420
720592	Serie_10908.txt	13167268	4a3ae7c490d7818f64199b55e298de94	doi:10.57745/VJ6DS2	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720592
721097	Serie_10909.txt	13205238	c0929e7568fed8dd5d790ef8ff5f17df	doi:10.57745/NJPIYM	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721097
720906	Serie_10910.txt	13215811	96d6660649de9a9e46e9f6be2c62605e	doi:10.57745/4ZD31V	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720906
721336	Serie_10911.txt	13196750	6b94dcc9f5164d09d34aec95021495ff	doi:10.57745/ABZYQO	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721336
720421	Serie_10912.txt	13194665	6d9791a927d60cbec8dacc3a8673433a	doi:10.57745/JRQGEG	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720421
720979	Serie_10913.txt	13198492	27c3494104b9fdd6217b47b0ba6a5c14	doi:10.57745/MH6YNW	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720979
720771	Serie_10914.txt	13191788	0ecf1e3469abcf728e1d2c3cc74ccc5e	doi:10.57745/EYKNH4	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720771
721096	Serie_10915.txt	13194595	9e03fba67ec4e2e814a0e8b1849b02a2	doi:10.57745/YPIJZB	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721096
721100	Serie_10916.txt	13205066	8eb28fdc4ff775874a2f883c183ec94a	doi:10.57745/K5IZFV	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721100
720924	Serie_10917.txt	13191123	d75818ff15410587bb9a884c84e0b726	doi:10.57745/SUO3SE	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720924
721329	Serie_10918.txt	13217283	4b5a4a3e42c68de608c06668e3ff6087	doi:10.57745/UFJ1SM	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721329
721259	Serie_10919.txt	13216164	d2eb1a3019577f77181cf1b6d148e88c	doi:10.57745/LXVPHX	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721259
720988	Serie_10920.txt	13206075	1563055e793ada4ee6e15ada59f4abdf	doi:10.57745/V6V3BW	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720988
721211	Serie_10921.txt	13209065	e7ac2379e2118221a16bd17dc8b5211f	doi:10.57745/8CZYYD	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721211
720682	Serie_10922.txt	13196659	e2fb352ffc886fdc999f648d5c82fc12	doi:10.57745/5IIR2Y	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720682
720471	Serie_10923.txt	13195393	d21ef926349486dabdf76c7695b57b9e	doi:10.57745/UBT2XM	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720471
720412	Serie_10924.txt	13215960	18fdc0398a0286a2e464f1a71e1ed23d	doi:10.57745/7HMD7B	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720412
721162	Serie_10925.txt	13203529	597d10aac2e7f69fdd2c5eae1972f19c	doi:10.57745/T50FKY	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721162
720970	Serie_10926.txt	13201843	faf58fe1925718827ed824bdf8d27291	doi:10.57745/SJQTQU	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720970
721158	Serie_10927.txt	13228077	a200d32be08b049f4389e1c9ff70160e	doi:10.57745/IVLGM3	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721158
721252	Serie_10928.txt	13196270	a58f930c1fb147192331acc58680526a	doi:10.57745/8PWOMK	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721252
720779	Serie_10929.txt	13215819	fc3a33f885d805b8f5b297db9a03be15	doi:10.57745/9FLNF9	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720779
721310	Serie_10930.txt	13231278	ccd3608cc4afb370e0210afd0642dc6b	doi:10.57745/BEYWUY	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721310
721196	Serie_10931.txt	13221520	f87cdd8052a2df9fef3086ce0ab10d5d	doi:10.57745/HRJROL	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721196
721345	Serie_10932.txt	13186816	f71128fdc892415cebccdf60e8ddc716	doi:10.57745/WODGXP	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721345
720936	Serie_10933.txt	13197201	d04f68af1181214f44d1e388049c286d	doi:10.57745/MDK0AR	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720936
720846	Serie_10934.txt	13200881	34dbee0af071cc90b471f32f77efede9	doi:10.57745/ZZV6FM	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720846
720788	Serie_10935.txt	13167123	f94142020f2c781c2748302c57fe5f89	doi:10.57745/VE1KSK	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720788
720792	Serie_10936.txt	13207412	216783066c2a65487bdacbe72e68de78	doi:10.57745/PQYE4E	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720792
721261	Serie_10937.txt	13202711	e08fcd3609fdf707c5c25f88ca5c15ec	doi:10.57745/KFCQD8	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721261
721063	Serie_10938.txt	13205907	c55f67b4aeb1397679d52d63aea54f8a	doi:10.57745/IZBOEC	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721063
720808	Serie_10939.txt	13227850	c23e074827494c7ee02a6d70c5e8beae	doi:10.57745/AUGLZK	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720808
720560	Serie_10940.txt	13199399	822d08fba905505c8daf7927b2c2152e	doi:10.57745/2TWMDD	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720560
720872	Serie_10941.txt	13208157	4a808590a4aa776f9cba5a48922a102b	doi:10.57745/TM1DEL	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720872
720966	Serie_10942.txt	13214783	0053af34c0781cda12be05cfa6e93832	doi:10.57745/H5DLVN	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720966
721248	Serie_10943.txt	13209529	a6daeb31d64e5252fe502ed4b53e4e98	doi:10.57745/II5IZ1	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721248
721081	Serie_10944.txt	13168034	f56aea592db5f5b6e62d11bd8a29c123	doi:10.57745/X2VZFG	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721081
720644	Serie_10945.txt	13196268	5c0808689ac22e3480f4beef8c85fdee	doi:10.57745/7PNDAQ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720644
720517	Serie_10946.txt	13218665	87091be9ca56c7727a0e7738b9f01f73	doi:10.57745/HCFFZG	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720517
721009	Serie_10947.txt	13212757	78ea206a75aa8a967821ce10af663116	doi:10.57745/JIDO3I	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721009
720845	Serie_10948.txt	13198115	a731cde9770d02581b2c4b2cb535fbbc	doi:10.57745/EZX5ZG	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720845
720689	Serie_10949.txt	13207451	ba02821e50c644f74ca34238bf2e0f53	doi:10.57745/RIDEF2	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720689
720647	Serie_10950.txt	13222992	99a376a6fd9490f39e72f7f100c3e60d	doi:10.57745/BDH4L9	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720647
720457	Serie_10951.txt	13200169	9c1e5b3c49e91791a0268a52215b82af	doi:10.57745/85ETTI	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720457
720892	Serie_10952.txt	13215478	dd23cf09e415fee9762cb2239ab574d3	doi:10.57745/RJAJSQ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720892
720939	Serie_10953.txt	13199718	dbbcb20cc43d966af19e367ead7f6f4d	doi:10.57745/7JVA6H	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720939
720907	Serie_10954.txt	13203987	5e8c32685e2c6e5290359c065752bf3d	doi:10.57745/4WQ6N0	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720907
721101	Serie_10955.txt	13198669	100eef483a2feb08b6670e0239ff4893	doi:10.57745/QUKR29	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721101
721138	Serie_10956.txt	13186712	ad4552c10ce0450459696592bc0843bc	doi:10.57745/VELNYU	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721138
720799	Serie_10957.txt	13185908	474f4fa44e4ae32b8f5a6591296e7d0a	doi:10.57745/XL6OUU	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720799
720509	Serie_10958.txt	13226174	897c248e87650eaeb439ccc3b8a37ac9	doi:10.57745/SUI82X	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720509
720669	Serie_10959.txt	13195635	932102499aa201c5a52071001842b0e2	doi:10.57745/NMDB2E	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720669
721028	Serie_10960.txt	13185628	44ab39a35dd173c0ba0a819fc3a4b2dd	doi:10.57745/NFDUMA	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721028
721012	Serie_10961.txt	13182628	8ac39ff003c8e6eef74e1cab0502a744	doi:10.57745/K9EAWO	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721012
720911	Serie_10962.txt	13204290	d4b0a77507a507a9f7256f098ea8d9bf	doi:10.57745/XN34XH	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720911
720417	Serie_10963.txt	13199616	893f581795063cfc58d84a0be7b2cf06	doi:10.57745/NXTHUJ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720417
720741	Serie_10964.txt	13219037	ad4416f55e47297c9d84efa5d0c98cd0	doi:10.57745/XKDT5U	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720741
720880	Serie_10965.txt	13218658	37b2c5440dddb8836104a6d1c6856dfb	doi:10.57745/FD3OCA	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720880
720931	Serie_10966.txt	13209121	caa315d4142f8d4e4c4a7bfe400bf013	doi:10.57745/JUFPPP	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720931
720640	Serie_10967.txt	13197827	6783bb7a403c908fc69e337d6fbf5430	doi:10.57745/SNGD9R	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720640
720409	Serie_10968.txt	13187562	7611a33f00d1b1ee0810c3280ece90f9	doi:10.57745/THBQZD	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720409
720540	Serie_10969.txt	13206873	edc87b708d9172ba6d005aa73e2579a3	doi:10.57745/H2KJP8	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720540
720650	Serie_10970.txt	13185652	a2c0547d55bf1d4ca1fcd010877210d5	doi:10.57745/OSV2Z4	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720650
720605	Serie_10971.txt	13188589	8459f041f1e31d21da493ec3ddb36476	doi:10.57745/DJHQJG	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720605
720764	Serie_10972.txt	13197119	cf4f5789f8caf52647d7884eb8f00968	doi:10.57745/YORKDQ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720764
720798	Serie_10973.txt	13229815	c59d405118fa91989be6db6c4ae21ff0	doi:10.57745/UONZDX	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720798
720405	Serie_10974.txt	13189222	12d7841d2df910a66fef7828e5d4455c	doi:10.57745/BGMWRA	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720405
721209	Serie_10975.txt	13208720	1f25860ca48f8075a8854f57508ff0d1	doi:10.57745/Z1RRWR	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721209
720775	Serie_10976.txt	13193660	a4dd6bdef155de61b4b0101dac60db4c	doi:10.57745/YO1OIL	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720775
720755	Serie_10977.txt	13218519	7289b1f09ab5188cc08ff1bed9b9c165	doi:10.57745/00LAZG	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720755
720484	Serie_10978.txt	13205915	6e8ae7fb5aba56a37eca3ad17977576c	doi:10.57745/ZYYNBW	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720484
720868	Serie_10979.txt	13210300	0e3d868b204db9999cc83423f4c2fb5f	doi:10.57745/SVDSPH	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720868
720710	Serie_10980.txt	13222612	3ac34eff6a2034aba403f20db1076dac	doi:10.57745/I8UVJR	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720710
720520	Serie_10981.txt	13237748	e87293354ce9868348aaca94b52253c4	doi:10.57745/Q6K1DM	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720520
721380	Serie_10982.txt	13208537	b8b84cdb19f68e538a71e632bb4613f1	doi:10.57745/LWJWCK	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721380
721371	Serie_10983.txt	13205094	40b831f8f899b35f96f08fd299ba99bd	doi:10.57745/JKLLNQ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721371
720863	Serie_10984.txt	13211417	e543b0b2bb37cfed061da7f6b339ec7b	doi:10.57745/EXYBAT	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720863
720504	Serie_10985.txt	13195695	c96f66d0ede4ca02ea88e13b2e53be20	doi:10.57745/M7TNER	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720504
721114	Serie_10986.txt	13201649	4bc2853f089b43c0d5cdff452647a943	doi:10.57745/BDBEDJ	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721114
720406	Serie_10987.txt	13193751	b737ef7b1156b6baf149d8626bbb1e02	doi:10.57745/SQXOYY	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720406
721316	Serie_10988.txt	13203459	409bd0d7e3da2cda9cce9baa17a1dc91	doi:10.57745/RCQ7ED	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721316
720628	Serie_10989.txt	13238232	b84ef8ee2eb5e04fb5c5e7807190465c	doi:10.57745/6K5KLV	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720628
720544	Serie_10990.txt	13190336	4aa0331de663a88fd2a3ba9d58bb671d	doi:10.57745/IBCSFE	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720544
720700	Serie_10991.txt	13185111	020216df5d87458db8d09aab3d8dda89	doi:10.57745/IGUOD7	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720700
721346	Serie_10992.txt	13223801	e5986af098a8e0c81bb8dc5de8fea297	doi:10.57745/WXS3LE	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721346
720598	Serie_10993.txt	13196075	d6f0255c9a13c492b9dd9a8b2bec8102	doi:10.57745/V5I9YY	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720598
721277	Serie_10994.txt	13202343	91324bc27f1fe9c5941d1f68a9daee23	doi:10.57745/FIEYXF	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721277
721082	Serie_10995.txt	13192588	f9498680e273844146730a4bc91d7e52	doi:10.57745/JXXJIS	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721082
720908	Serie_10996.txt	13203882	f2a3f1e8ae5bdf3466340d56d7b258c1	doi:10.57745/QK8UKU	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720908
720900	Serie_10997.txt	13172162	c98ff466f459d33cf09f6239e0390e52	doi:10.57745/X6MNK0	https://entrepot.recherche.data.gouv.fr/api/access/datafile/720900
721029	Serie_10998.txt	13208427	418446757a16b4f2044145694b7eea0e	doi:10.57745/P3QLXY	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721029
721086	Serie_10999.txt	13216581	34cfb2ff14efdf2cd049a2ec4fc766c5	doi:10.57745/0KEKUM	https://entrepot.recherche.data.gouv.fr/api/access/datafile/721086
~~~
