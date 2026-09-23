import {chromium,expect} from '@playwright/test';
import {mkdir} from 'node:fs/promises';
const base=process.env.HUNDUN_UI_URL||'http://127.0.0.1:8766';
const browser=await chromium.launch({headless:true,...(process.env.HUNDUN_BROWSER?{executablePath:process.env.HUNDUN_BROWSER}:{})});
const page=await browser.newPage({viewport:{width:1440,height:1000}});const errors=[];page.on('pageerror',e=>errors.push(e.message));
await mkdir('check/browser',{recursive:true});
try{
 await page.goto(base);await expect(page.getByLabel('计算主机')).toBeVisible();
 await page.getByLabel('连接与模型').click();await expect(page.getByRole('dialog')).toBeVisible();
 await page.getByLabel('名称',{exact:true}).fill('本机测试');await page.getByRole('button',{name:'保存连接',exact:true}).click();
 await expect(page.getByRole('status')).toContainText('已保存');await page.getByRole('button',{name:'探测能力'}).click();await expect(page.getByRole('status')).toContainText('连接已就绪');
 await page.screenshot({path:'check/browser/settings.png'});
 await page.getByLabel('关闭面板').click();await page.getByRole('button',{name:'AI 助手',exact:true}).click();
 await page.getByLabel('模拟目标或追问').fill('检查可用工作目录');await page.getByRole('button',{name:'开始任务',exact:true}).click();
 await expect(page.locator('.client-taskinfo')).toContainText('配置模型后继续',{timeout:15000});
 await expect(page.locator('.client-events')).toContainText('配置模型端点');
 await page.getByRole('button',{name:'取消 AI 任务',exact:true}).click();await expect(page.locator('.client-taskinfo')).toContainText('已取消',{timeout:10000});
 await page.screenshot({path:'check/browser/agent.png'});
 await page.getByLabel('关闭面板').click();await page.getByRole('button',{name:'输入版本',exact:true}).click();await expect(page.getByLabel('导入目录')).toBeVisible();await page.keyboard.press('Escape');
 for(const name of ['算例','网格与边界','计算监看','后处理','报告中心']){await page.getByRole('navigation').getByRole('button',{name,exact:true}).click();await expect(page.getByRole('heading',{name,exact:true,level:1})).toBeVisible()}
 await page.screenshot({path:'check/browser/desktop.png',fullPage:true});
 await page.setViewportSize({width:390,height:844});await page.getByRole('button',{name:'AI 助手',exact:true}).click();await expect(page.getByRole('dialog')).toBeVisible();
 expect(await page.evaluate(()=>document.documentElement.scrollWidth<=window.innerWidth)).toBeTruthy();await page.screenshot({path:'check/browser/mobile.png'});
 expect(errors).toEqual([]);console.log('PASS: settings, host probe, task SSE, cancel, inputs, six pages, mobile; console errors=0');
}finally{await browser.close()}
