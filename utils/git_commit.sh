#! /bin/sh
# run this script at the root of a git repo to extract its latest git commit

if [ "$1" = "--debug" ]; then
    _debug_mode=1
fi

GIT_STATE="unknown"
GITHASH="?"

if command -v git >/dev/null 2>&1; then
    [ -z "$_debug_mode" ] || echo "git detected!"
    _repo_status=$(git branch -r 1>/dev/null 2>&1 && echo "0")

    if [ "$_repo_status" = "0" ]; then
                
        GIT_STATUS_LINES=$(git status --short -uno | wc -l)
        if [ "$GIT_STATUS_LINES" = 0 ]; then
            GIT_STATE="unmodified"
        else
            GIT_STATE="modified"
        fi
        GITHASH=$(git rev-parse HEAD)

    fi
else
    [ -z "$_debug_mode" ] || echo "extracting information without git!"
    if [ -f .git/HEAD ]; then
        _head_cat=$(cat .git/HEAD)
        BRANCH=${_head_cat#ref: }
        if [ -f ".git/$BRANCH" ]; then
            GITHASH=$(cat ".git/$BRANCH")
        fi
    fi
fi

# last resort: try extracting git hash from environment
if [ "$GITHASH" = "?" ]; then
    [ -z "$_debug_mode" ] || echo "extracting git reposiroty information from env!"
    if [ -n "$CI_COMMIT_SHA" ]; then # gitlab ci/woodpecker ci
       GITHASH="$CI_COMMIT_SHA (from CI_COMMIT_SHA)"
    fi
fi


echo "$GITHASH (repo status: $GIT_STATE)"