# Fixture plan

## 3. Phases

### Phase 1: First

| Step | What | Agent | Lane | Depends | Origin |
|---|---|---|---|---|---|
| 1.1 | A grandfathered step. | orch | | | |
| 4.1 | A step that closed. | worker | tools | | |

### Phase 4: Second

| Step | What | Agent | Lane | Depends | Origin |
|---|---|---|---|---|---|
| 4.2 | A step that closes in this commit. | worker | tools | | |
| 4.3 | An open step. | worker | tools | | |

## 4. Status

| Step | Status | Commit | Duration |
|---|---|---|---|
| 1.1 | done: grandfathered, no duration | this plan's commit | |
| 4.1 | done: closed; world-0088 | aaaaaaa, bbbbbbb | worker 12 min; 40 min from ccccccc |
| 4.2 | part: closes here | this commit | 1 h 5 min from bbbbbbb |

## 5. For the user

- Nothing.
