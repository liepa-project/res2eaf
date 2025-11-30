# Orginal script
```
uv run res2eaf.py -j -l test/test.lat -o target/test_org.eaf
```

# Packaged script
```
uv run res2eaf -j -l test/test.lat -o target/test_lib.eaf
diff target/test_lib.eaf target/test_org.eaf
```