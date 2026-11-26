"""从 replay 下已登记的文件重建公开复算包，各成员使用固定时间戳。"""
from pathlib import Path
import subprocess
from tempfile import NamedTemporaryFile
from zipfile import ZIP_DEFLATED, ZipFile, ZipInfo

root = Path(__file__).resolve().parent
source = root / "replay"
paths = subprocess.check_output(["git", "ls-files", "-z", "--", "."], cwd=source)
names = sorted(item.decode("utf-8") for item in paths.split(b"\0") if item)
if not names:
    raise RuntimeError("replay 目录缺少已登记的源码")
archive = root / "ch39-omega-math-replay.zip"
with NamedTemporaryFile(dir=root, suffix=".zip", delete=False) as tmp:
    candidate = Path(tmp.name)
try:
    with ZipFile(candidate, "w") as output:
        for name in names:
            member = ZipInfo(name)
            member.create_system = 3
            member.external_attr = 0o100644 << 16
            output.writestr(member, (source / name).read_bytes(), compress_type=ZIP_DEFLATED)
    with ZipFile(candidate) as output:
        if output.namelist() != names or any(
            output.read(name) != (source / name).read_bytes() for name in names
        ):
            raise RuntimeError("归档内容与复算目录不一致")
    candidate.chmod(0o644)
    candidate.replace(archive)
finally:
    candidate.unlink(missing_ok=True)
print(f"已核对 {len(names)} 个文件")
