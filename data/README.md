### Example of generating completely random data of length 10000 with destination test.txt
```python generate_data.py --random-type=normal --filename=test.txt --seq-len=10000```

### Example of random seq of length 500 with pattern of length 10 repeating 20 times
```python generate_data.py --random-type=repeating --filename=test_repeating.txt --seq-len=500 --pattern-len=10 --pattern-freq=20```

### Given file input.txt, invoke edit mode = true and give number of inserts/deletes/swaps, write to file with format input_edited_<x>I_<y>D_<z>S.txt where x,y,z is num ins/dels/swps respectively
```python generate_data.py --filename=input.txt --edit-mode=True --num-inserts=4 --num-deletes=6 --num-swaps=3```
