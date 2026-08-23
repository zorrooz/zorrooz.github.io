/**
 * Markdown 渲染后的 DOM 增强：标题锚点、代码块包装（语言标签 + 复制）、表格复制。
 * 全部为幂等操作（已处理元素跳过），在每次内容渲染后执行。
 */
import i18n from '@/i18n/index'
import { iconSvg } from '@/utils/icons'
import { scrollToHeading } from '@/utils/scroll'
import { HEADER_OFFSET } from '@/config'

const t = (key: string) => i18n.global.t(key)

/** 复制回调：由组件注入（带按钮反馈） */
export type CopyHandler = (text: string, button: HTMLButtonElement) => void

/** 标题增强：清理与页面标题重复的 h1，并为全部标题挂 # 锚点按钮 */
export function enhanceHeadings(container: HTMLElement, articleTitle: string) {
  if (articleTitle) {
    const pageTitle = articleTitle.trim().toLowerCase()
    container.querySelectorAll('h1').forEach((h1) => {
      const h1Text = h1.textContent.trim().toLowerCase()
      if (h1Text === pageTitle) {
        h1.remove()
        return
      }
      const h2 = document.createElement('h2')
      Array.from(h1.attributes).forEach((attr) => h2.setAttribute(attr.name, attr.value))
      h2.innerHTML = h1.innerHTML
      h1.replaceWith(h2)
    })
  }

  container.querySelectorAll('h2, h3, h4, h5, h6').forEach((heading) => {
    heading.querySelector('.heading-anchor')?.remove()
    const anchorBtn = Object.assign(document.createElement('button'), {
      type: 'button',
      className: 'heading-anchor',
      textContent: '#',
      ariaLabel: t('anchorHeading'),
      tabIndex: 0,
      ariaHidden: 'false',
    })

    const jump = () => {
      scrollToHeading(heading, HEADER_OFFSET)
      setTimeout(() => anchorBtn.blur(), 300)
    }
    anchorBtn.addEventListener('click', (e) => {
      e.stopPropagation()
      jump()
    })
    anchorBtn.addEventListener('keydown', (e) => {
      if (e.key === 'Enter' || e.key === ' ') {
        e.preventDefault()
        jump()
      }
    })

    heading.appendChild(anchorBtn)
  })
}

/** 代码块增强：包 wrapper + 语言标签头 + 复制按钮 */
export function enhanceCodeBlocks(container: HTMLElement, onCopy: CopyHandler) {
  container.querySelectorAll('pre').forEach((pre) => {
    if (pre.querySelector('.code-block-header')) return
    const code = pre.querySelector('code')
    if (!code) return

    const language = (code.className.match(/language-(\w+)/) || ['', 'text'])[1]
    const header = document.createElement('div')
    header.className = 'code-block-header d-flex align-items-center justify-content-between'

    const langLabel = document.createElement('span')
    langLabel.className = 'code-language'
    langLabel.textContent = language

    const copyButton = document.createElement('button')
    copyButton.type = 'button'
    copyButton.className = 'copy-button btn-icon d-flex align-items-center justify-content-center'
    copyButton.setAttribute('aria-label', t('copyCode'))
    copyButton.innerHTML = iconSvg('copy', 16)
    copyButton.addEventListener('click', () => onCopy(code.textContent ?? '', copyButton))

    header.append(langLabel, copyButton)
    const wrapper = document.createElement('div')
    wrapper.className = 'code-block-wrapper'
    pre.parentNode?.insertBefore(wrapper, pre)
    wrapper.append(header, pre)
  })
}

/** 表格增强：包 wrapper + 悬停复制按钮（TSV 格式，可直接粘贴 Excel） */
export function enhanceTables(container: HTMLElement, onCopy: CopyHandler) {
  container.querySelectorAll('table').forEach((table) => {
    if (table.closest('.table-copyable')) return

    const copyButton = document.createElement('button')
    copyButton.type = 'button'
    copyButton.className = 'table-copy-btn'
    copyButton.setAttribute('aria-label', t('copyTable'))
    copyButton.innerHTML = iconSvg('copy', 16)
    copyButton.addEventListener('click', () => onCopy(tableToTsv(table), copyButton))

    const wrapper = document.createElement('div')
    wrapper.className = 'table-copyable'
    table.parentNode?.insertBefore(wrapper, table)
    wrapper.append(copyButton)
    wrapper.append(table)
  })
}

function tableToTsv(table: HTMLTableElement): string {
  return Array.from(table.querySelectorAll('tr'))
    .map((tr) =>
      Array.from(tr.querySelectorAll('th, td'))
        .map((cell) => (cell.textContent || '').trim().replace(/\s+/g, ' '))
        .join('\t'),
    )
    .join('\n')
}
