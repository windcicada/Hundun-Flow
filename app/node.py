"""Versioned SSH helper. JSON on stdin/stdout; native processes stay detached."""

import json
import sys
from .ops import Operations
from .jobs import ApiError
from .data import safe


def main():
    try:
        p = json.load(sys.stdin)
        result = Operations(p["profile"], p["state"]).dispatch(
            p["action"], p.get("params", {})
        )
        print(json.dumps({"result": safe(result)}))
    except ApiError as exc:
        print(
            json.dumps(
                {
                    "error": {
                        "code": exc.code,
                        "message": exc.message,
                        "status": exc.status,
                    }
                }
            )
        )
    except Exception as exc:
        print(
            json.dumps(
                {"error": {"code": "node_failed", "message": str(exc), "status": 409}}
            )
        )


if __name__ == "__main__":
    main()
