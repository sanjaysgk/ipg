#!/usr/bin/env bats

setup() {
    REPO_ROOT="$(git rev-parse --show-toplevel)"
    SCRIPT="${REPO_ROOT}/bin/fetch_test_bundle.sh"
    TMPDIR_REAL="$(mktemp -d)"
    export IPG_TEST_BUNDLE_DIR="${TMPDIR_REAL}/.test-bundle"
}

teardown() {
    rm -rf "${TMPDIR_REAL}"
}

@test "script exists and is executable" {
    [ -x "${SCRIPT}" ]
}

@test "fails with clear error on bad URL" {
    run env IPG_TEST_BUNDLE_URL="https://example.invalid/nope.tar.gz" "${SCRIPT}"
    [ "${status}" -ne 0 ]
    [[ "${output}" =~ "download failed" ]]
}

@test "skips download if .test-bundle already populated" {
    mkdir -p "${IPG_TEST_BUNDLE_DIR}/reference"
    touch "${IPG_TEST_BUNDLE_DIR}/samplesheet_test.csv"
    run "${SCRIPT}"
    [ "${status}" -eq 0 ]
    [[ "${output}" =~ "already present" ]]
}
