"""Installed client entry point."""

import argparse
from platformdirs import user_state_dir
from . import __version__


def main():
    parser = argparse.ArgumentParser(prog="hundun-client")
    parser.add_argument("--version", action="version", version=__version__)
    commands = parser.add_subparsers(dest="command", required=True)
    serve = commands.add_parser("serve", help="启动本机工作台")
    serve.add_argument("--port", type=int, default=8765)
    serve.add_argument("--state", default=user_state_dir("hundun-client"))
    args = parser.parse_args()
    if args.command == "serve":
        import uvicorn
        from .backend import create_app

        uvicorn.run(create_app(state=args.state), host="127.0.0.1", port=args.port)


if __name__ == "__main__":
    main()
