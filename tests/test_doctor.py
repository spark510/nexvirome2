import subprocess
from nexvirome2 import doctor
from nexvirome2.cli import main


def test_missing_tools_exit_nonzero(monkeypatch, capsys):
    monkeypatch.setattr(doctor.shutil, 'which', lambda _: None)
    assert main(['doctor']) == 1
    assert 'not found in PATH' in capsys.readouterr().out


def test_present_but_broken_tool_not_healthy(monkeypatch):
    monkeypatch.setattr(doctor.shutil, 'which', lambda name: '/env/bin/'+name)
    monkeypatch.setattr(doctor.subprocess, 'run', lambda args, **kwargs:
                        subprocess.CompletedProcess(args, 1, '', 'missing shared library'))
    result = doctor.check()
    assert not result['ok']
    assert all(row['returncode'] == 1 for row in result['tools'])


def test_all_tools_start(monkeypatch):
    monkeypatch.setattr(doctor.shutil, 'which', lambda name: '/env/bin/'+name)
    monkeypatch.setattr(doctor.subprocess, 'run', lambda args, **kwargs:
                        subprocess.CompletedProcess(args, 0, 'test version', ''))
    assert doctor.check()['ok']
