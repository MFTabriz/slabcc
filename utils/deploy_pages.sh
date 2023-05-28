#! /bin/bash
# shellcheck disable=SC2086

# by default, this script commits all the files in "public" directory to "pages" branch of the present repository keeping only latest 10 commits
# all commandline flags get passed to "cp" and "git clone" commands. typically it is set to "--verbose" for testing
# options:
#
## $PAGES_PASS (required): authorization token for pushing commits to the origin repository
## $STAGE_PATH           : path to the directory which contains files to be deployed (default: 'public')
## $HISTORY_LENGTH       : max number of latest commits to keep (default: 10)
## $PAGES_USER           : username for authorization (default: original commit author)


# todo:
#       make a plugin for woodpecker-ci
#       for tag events, tag also the commit in pages branch


set -e

color_red="\033[31;49;1m"
color_green="\033[32;49;1m"
color_yellow="\033[33;49;1m"
color_reset="\033[0m"


if [ -z "$STAGE_PATH" ]; then 
    STAGE_PATH="public"
fi

if ! [ -d "$STAGE_PATH" ]; then 
    printf "%bERROR: CANNOT ACCESS STAGE DIRECTORY: $STAGE_PATH%b\n" "${color_red}" "${color_reset}"
    exit 1
fi

if [ -z "$PAGES_PASS" ]; then
    printf "%bWARNING: \$PAGES_PASS is not set! Skipping deployment to pages%b\n" "${color_yellow}" "${color_reset}"
    tar -czvf pages.tar.gz "$STAGE_PATH"/
    curl -s -v --upload-file pages.tar.gz https://transfer.sh/pages.tar.gz 2>&1 | tee .upload.log
    full_link=$(grep "x-url-delete:" .upload.log | sed 's/x-url-delete://')
    download_link=${full_link/%pages.tar.gz*/pages.tar.gz}
    echo -e "\n\nDOWNLOAD LINK:$download_link\n"
    sha256sum pages.tar.gz
    exit 0
fi

if [ -z "$HISTORY_LENGTH" ]; then 
    HISTORY_LENGTH=10
fi

if [ -z "$PAGES_USER" ]; then 
    PAGES_USER="${CI_COMMIT_AUTHOR}"
fi

if [ -n "$1" ]; then # pass extra flags and activate verbose mode
    EXTRA_FLAGS=( "$@" )
    VERBOSE_COMMIT="--allow-empty"
    set -x
fi

DEPLOY_BRANCH="pages"
STAGE_FULL_PATH=$(realpath "$STAGE_PATH")

echo "Staging: $STAGE_FULL_PATH" && ls -lah "$STAGE_PATH"
echo "Max history length: $HISTORY_LENGTH"
echo "Extra flags: ${EXTRA_FLAGS[*]}"
MAX_OLD_COMMITS=$(("$HISTORY_LENGTH" - 1))

TMP_DIR=$(mktemp -d)
TMP_BRANCH=$(tr -dc A-Za-z0-9 < /dev/urandom | head -c 20)
# clone the DEPLOY_BRANCH if available; otherwise, just clone the head of main branch

git clone "${EXTRA_FLAGS[@]}" --single-branch --branch="$DEPLOY_BRANCH" --depth="$MAX_OLD_COMMITS" -- https://"${PAGES_USER}":"${PAGES_PASS}"@"${CI_REPO_REMOTE##https://}" "$TMP_DIR" \
|| git clone "${EXTRA_FLAGS[@]}" --depth=1 -- https://"${PAGES_USER}":"${PAGES_PASS}"@"${CI_REPO_REMOTE##https://}" "$TMP_DIR"

cd "$TMP_DIR" 

git config user.email "${CI_COMMIT_AUTHOR}@noreply.codeberg.org" && git config user.name "${CI_COMMIT_AUTHOR}"

# cleanup the git repository, keep only the latest commits if DEPLOY_BRANCH already exists
( 
    git switch "${DEPLOY_BRANCH}" && NEW_GIT_BASE=$(git rev-list --max-parents=0 HEAD) &&\
    git checkout --orphan "${TMP_BRANCH}" "$NEW_GIT_BASE" &&\
    git commit --reuse-message="$NEW_GIT_BASE" &&\
    git rebase --onto "${TMP_BRANCH}" "$NEW_GIT_BASE" "$DEPLOY_BRANCH" &&\
    git branch -D "${TMP_BRANCH}" &&\
    git rm -fr . && \
    git checkout HEAD -- .woodpecker.yml &&\ # workaround https://codeberg.org/Codeberg-CI/feedback/issues/114
    git reflog expire --expire=all --all && git prune --progress && git gc --aggressive
) || git switch --orphan "$DEPLOY_BRANCH"

# add new files and commit
cp -r "${EXTRA_FLAGS[@]}" "${STAGE_FULL_PATH}"/* ./
git add .

if ! git commit -m "deploying $DEPLOY_BRANCH" -m "source commit $CI_COMMIT_SHA" -m "Co-authored-by: deploy-bot <deploy-bot@noreply.codeberg.org>" ${VERBOSE_COMMIT}; then
    printf "%bWARNING: git commit failed! probably there's nothing new to deploy.%b\n" "${color_yellow}" "${color_reset}" 
    exit 0
fi

# push
git push -f || git push --set-upstream origin "$DEPLOY_BRANCH"
git log --decorate --oneline

printf "%bDONE: %s was successfully deployed to %s branch!%b\n" "${color_green}" "$STAGE_PATH" "${DEPLOY_BRANCH}" "${color_reset}"
echo "Pages should be accessible at: https://${CI_REPO_OWNER}.codeberg.page/${CI_REPO_NAME}"