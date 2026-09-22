// SPDX-License-Identifier: Apache-2.0
// Read-only browser smoke: run against an already started local workbench.
import { chromium, expect } from '@playwright/test';
import { mkdir } from 'node:fs/promises';
import { fileURLToPath } from 'node:url';

const base = process.env.HUNDUN_UI_URL || 'http://127.0.0.1:8765';
const output = fileURLToPath(new URL('./.qa/', import.meta.url));
await mkdir(output, {recursive: true});
const browser = await chromium.launch({
  headless: true,
  ...(process.env.HUNDUN_BROWSER ? {executablePath: process.env.HUNDUN_BROWSER} : {}),
});
const page = await browser.newPage({viewport: {width: 1600, height: 1060}});
const errors = [];
page.on('pageerror', error => errors.push(error.message));
try {
  await page.goto(base);
  await expect(page.getByLabel('当前运行')).toBeVisible();
  const response = await page.request.get(base + '/api/runs');
  expect(response.ok()).toBeTruthy();
  const {runs} = await response.json();
  const run = runs.find(value => value.seconds_per_step != null);
  if (run) {
    await page.getByLabel('当前运行').selectOption(run.id);
    await expect(page.locator('.metric-value').first()).not.toBeEmpty();
  }
  await page.screenshot({path: output + 'monitor.png', fullPage: true});
  for (const name of ['算例', '网格与边界', '计算监看', '后处理', '报告中心']) {
    await page.getByRole('navigation').getByRole('button', {name, exact: true}).click();
    await expect(page.getByRole('heading', {level: 1, name, exact: true})).toBeVisible();
  }
  if (run) {
    const report = page.waitForEvent('download');
    await page.getByRole('button', {name: '导出报告', exact: true}).click();
    await (await report).saveAs(output + 'report.md');
  }
  // Field rendering is opt-in because a production frame may be large.
  if (process.env.HUNDUN_TEST_FIELD_RUN) {
    await page.getByRole('navigation').getByRole('button', {name: '后处理', exact: true}).click();
    await page.getByLabel('当前运行').selectOption(process.env.HUNDUN_TEST_FIELD_RUN);
    await expect(page.getByRole('button', {name: '生成截面', exact: true})).toBeEnabled({timeout: 30000});
    await page.getByRole('button', {name: '生成截面', exact: true}).click();
    await expect(page.locator('.field-canvas img')).toBeVisible({timeout: 120000});
    await page.screenshot({path: output + 'field.png', fullPage: true});
  }
  await page.setViewportSize({width: 390, height: 844});
  await page.waitForTimeout(350);
  expect(await page.evaluate(() => document.documentElement.scrollWidth <= innerWidth)).toBeTruthy();
  await page.getByRole('button', {name: '打开导航', exact: true}).click();
  await page.getByRole('navigation').getByRole('button', {name: '工作台', exact: true}).click();
  await expect(page.locator('.sidebar')).not.toHaveClass(/open/);
  await page.waitForTimeout(350);
  await page.screenshot({path: output + 'mobile.png', fullPage: true});
  expect(errors).toEqual([]);
  console.log(JSON.stringify({ok: true, runs: runs.length, field: process.env.HUNDUN_TEST_FIELD_RUN || null, screenshots: output}));
} finally {
  await browser.close();
}
