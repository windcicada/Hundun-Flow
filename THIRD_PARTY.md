# Third-party software

Hundun Client is licensed under Apache-2.0. Python dependencies are installed as separate distributions, each retaining its upstream metadata and license files. The exact tested environment is recorded in `docs/deps.tsv`; installation requirements are in `pyproject.toml`.

| Component | Use | License |
|---|---|---|
| nanobot-ai 0.3.5 | Agent SDK and session runtime | MIT |
| MCP Python SDK | Typed simulation tools | MIT |
| FastAPI, Uvicorn, HTTPX, Pydantic | Local API and validation | MIT / BSD-3-Clause |
| AsyncSSH | SSH and SFTP adapter | EPL-2.0 OR GPL-2.0-or-later |
| platformdirs, keyring | Platform state and credentials | MIT |
| NumPy, Matplotlib, PyVista, VTK | Field reading and scientific plotting | BSD / PSF-based |
| React, React DOM | Bundled browser application | MIT |
| Lucide | Bundled browser icons | ISC |
| TypeScript, Vite, Playwright | Development and browser checks | Apache-2.0 / MIT |

Bundled frontend attribution and license text are in `docs/licenses.txt`. Python wheels retain complete upstream licensing in their distribution metadata. Solver distributions selected by the user retain their own licenses.
