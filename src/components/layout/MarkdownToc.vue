<!-- MarkdownTpc.vue-->
<template>
  <div class="toc-item sidebar-item">
    <h4 class="toc-title">目录</h4>
    <ul class="toc-list">
      <li v-for="item in tocList" :key="item.id" :style="{ paddingLeft: (item.level - 2) * 15 + 'px' }">
        <a href="javascript:;" @click="scrollToHeading(item.id)" class="toc-link">
          {{ item.text }}
        </a>
      </li>
    </ul>
  </div>
</template>

<script>
export default {
  props: {
    containerSelector: {
      type: String,
      required: true,
      default: '.markdown-body'
    },
    levels: {
      type: Array,
      default: () => [2, 3, 4, 5, 6]
    }
  },
  data() {
    return { tocList: [] }
  },
  mounted() {
    this.generateToc()
  },
  methods: {
    generateToc() {
      const tryGetHeadings = () => {
        const contentEl = document.querySelector(this.containerSelector)
        if (!contentEl) {
          console.warn(`未找到文章容器: ${this.containerSelector}`)
          setTimeout(tryGetHeadings, 500)
          return
        }

        const selectors = this.levels.map(l => `h${l}`).join(', ')
        const headings = Array.from(contentEl.querySelectorAll(selectors))

        if (headings.length === 0) {
          setTimeout(tryGetHeadings, 800)
          return
        }

        this.tocList = headings.map(heading => {
          let id = heading.id
          if (!id) {
            id = heading.textContent.trim()
              .toLowerCase()
              .replace(/[^\w\u4e00-\u9fa5-]/g, '-')
              .replace(/-+/g, '-')
              .replace(/^-|-$/g, '')
            heading.id = id
          }
          return {
            id,
            text: heading.textContent.trim(),
            level: Number(heading.tagName.replace('H', ''))
          }
        })
      }

      tryGetHeadings()
    },

    scrollToHeading(id) {
      const target = document.getElementById(id)
      if (!target) return

      const headerHeight = document.querySelector('header')?.offsetHeight || 0
      const targetTop = target.getBoundingClientRect().top + window.scrollY - headerHeight - 10

      window.scrollTo({
        top: targetTop,
        behavior: 'smooth'
      })
    }
  }
}
</script>

<style scoped>
.toc-title {
  margin: 0 0 12px;
  font-size: 15px;
  color: #333;
}

.toc-list {
  list-style: none;
  margin: 0;
  padding: 0;
}

.toc-link {
  display: block;
  text-decoration: none;
  color: #409eff;
  font-size: 14px;
  line-height: 1.8;
  transition: color 0.2s ease;
  white-space: nowrap;
  overflow: hidden;
  text-overflow: ellipsis;
}

.toc-link:hover {
  color: #66b1ff;
}

.toc-empty {
  margin: 0;
  font-size: 14px;
  color: #999;
  line-height: 1.5;
}
</style>
