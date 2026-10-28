"""qn_util に付属するサンプルスクリプトを作業ディレクトリへ配置する。"""

import argparse
import shutil
from importlib.resources import files
from pathlib import Path


def install(force: bool = False, destination: Path | None = None) -> list[Path]:
    """同梱 scripts の全ファイルを指定先（既定はカレント）へコピーする。"""
    target = Path.cwd() if destination is None else Path(destination)
    source = files("qn_util").joinpath("scripts")
    if not source.is_dir():
        raise FileNotFoundError("パッケージ内に scripts がありません。配布設定を確認してください。")
    items = list(source.iterdir())
    conflicts = [target / item.name for item in items if (target / item.name).exists()]
    if conflicts and not force:
        names = ", ".join(str(x) for x in conflicts)
        raise FileExistsError(f"同名ファイルが存在します: {names}（上書きには --force を指定）")
    copied = []
    for item in items:
        dst = target / item.name
        if item.is_dir():
            shutil.copytree(str(item), dst, dirs_exist_ok=force)
        else:
            with item.open("rb") as src, dst.open("wb") as out:
                shutil.copyfileobj(src, out)
        copied.append(dst)
    return copied


def main() -> None:
    """qn_util コマンドのエントリーポイント。"""
    parser = argparse.ArgumentParser(prog="qn_util")
    sub = parser.add_subparsers(dest="command", required=True)
    cmd = sub.add_parser("install", help="付属スクリプトをカレントディレクトリへコピー")
    cmd.add_argument("--force", action="store_true", help="同名ファイルを上書き")
    args = parser.parse_args()
    if args.command == "install":
        try:
            copied = install(force=args.force)
        except (FileExistsError, FileNotFoundError) as exc:
            parser.exit(1, f"エラー: {exc}\n")
        for path in copied:
            print(f"コピー: {path}")


if __name__ == "__main__":
    main()
