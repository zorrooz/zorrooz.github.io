<template>
  <div class="article-tree">
    <!-- 二级分类标题 -->
    <h2 class="category-title">{{ currentCategory.title }}</h2>

    <!-- 分类描述（可选） -->
    <p v-if="currentCategory.desc" class="category-desc">{{ currentCategory.desc }}</p>

    <!-- 文章列表 -->
    <div class="article-list">
      <div v-for="item in currentCategory.items" :key="item.name" class="article-item"
        :class="{ 'external-link': item.type === 'external' }" @click="handleItemClick(item)">
        <div class="article-header">
          <h3 class="article-title">{{ item.name }}</h3>
          <span v-if="item.type === 'external'" class="external-badge" title="外部链接">
            <svg xmlns="http://www.w3.org/2000/svg" width="12" height="12" viewBox="0 0 24 24" fill="none"
              stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
              <path d="M18 13v6a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2V8a2 2 0 0 1 2-2h6"></path>
              <polyline points="15 3 21 3 21 9"></polyline>
              <line x1="10" y1="14" x2="21" y2="3"></line>
            </svg>
          </span>
        </div>

        <p class="article-desc">{{ item.desc }}</p>

        <div class="article-footer">
          <span class="article-type" :class="item.type">
            {{ item.type === 'internal' ? '内部资源' : '外部资源' }}
          </span>
          <button class="view-btn" :class="{ 'external-btn': item.type === 'external' }">
            {{ item.type === 'internal' ? '查看详情' : '访问链接' }}
          </button>
        </div>
      </div>
    </div>

    <!-- 更新提示 -->
    <div v-if="showUpdateTip" class="update-tip">
      <p>{{ categoryHelper.updateTip }}</p>
      <router-link v-if="categoryHelper.updateUrlType === 'internal'" :to="categoryHelper.updateUrl"
        class="update-link">
        查看更新
      </router-link>
      <a v-else :href="categoryHelper.updateUrl" target="_blank" rel="noopener noreferrer" class="update-link">
        查看更新
      </a>
    </div>
  </div>
</template>

<script>
export default {
  name: 'NavigationTree',
  props: {
    // 接收一个二级分类对象（来自categoryList中的一个对象）
    category: {
      type: Object,
      required: true,
      validator(value) {
        return value.title && Array.isArray(value.items);
      }
    },
    // 接收categoryHelper对象
    categoryHelper: {
      type: Object,
      required: true
    },
    // 是否显示更新提示
    showUpdateTip: {
      type: Boolean,
      default: true
    }
  },
  computed: {
    currentCategory() {
      return this.category;
    }
  },
  methods: {
    handleItemClick(item) {
      if (item.type === 'internal') {
        // 内部链接，使用Vue Router导航
        this.$router.push(item.url);
      } else {
        // 外部链接，在新窗口打开
        window.open(item.url, '_blank', 'noopener,noreferrer');
      }

      // 发出点击事件，供父组件监听
      this.$emit('item-click', item);
    }
  }
};
</script>

<style scoped>
.article-tree {
  max-width: 800px;
  margin: 0 auto;
  padding: 20px;
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
}

.category-title {
  font-size: 1.8rem;
  font-weight: 600;
  color: #2c3e50;
  margin-bottom: 12px;
  padding-bottom: 8px;
  border-bottom: 2px solid #3498db;
}

.category-desc {
  font-size: 1rem;
  color: #7f8c8d;
  margin-bottom: 24px;
  line-height: 1.6;
}

.article-list {
  display: flex;
  flex-direction: column;
  gap: 20px;
}

.article-item {
  background: #ffffff;
  border: 1px solid #e0e0e0;
  border-radius: 8px;
  padding: 20px;
  transition: all 0.3s ease;
  cursor: pointer;
}

.article-item:hover {
  transform: translateY(-2px);
  box-shadow: 0 4px 12px rgba(0, 0, 0, 0.1);
  border-color: #3498db;
}

.article-item.external-link {
  border-left: 4px solid #e74c3c;
}

.article-item.external-link:hover {
  border-color: #e74c3c;
}

.article-header {
  display: flex;
  align-items: center;
  justify-content: space-between;
  margin-bottom: 12px;
}

.article-title {
  font-size: 1.25rem;
  font-weight: 600;
  color: #2c3e50;
  margin: 0;
  line-height: 1.4;
}

.external-badge {
  display: inline-flex;
  align-items: center;
  justify-content: center;
  width: 24px;
  height: 24px;
  background: #e74c3c;
  color: white;
  border-radius: 50%;
  cursor: help;
}

.article-desc {
  font-size: 0.95rem;
  color: #555;
  margin: 16px 0;
  line-height: 1.6;
}

.article-footer {
  display: flex;
  align-items: center;
  justify-content: space-between;
  margin-top: 16px;
  padding-top: 16px;
  border-top: 1px solid #f0f0f0;
}

.article-type {
  font-size: 0.85rem;
  padding: 4px 8px;
  border-radius: 4px;
  font-weight: 500;
}

.article-type.internal {
  background: #d5e8f8;
  color: #2980b9;
}

.article-type.external {
  background: #fadbd8;
  color: #c0392b;
}

.view-btn {
  padding: 8px 16px;
  background: #3498db;
  color: white;
  border: none;
  border-radius: 4px;
  cursor: pointer;
  font-size: 0.9rem;
  font-weight: 500;
  transition: background-color 0.2s ease;
}

.view-btn:hover {
  background: #2980b9;
}

.view-btn.external-btn {
  background: #e74c3c;
}

.view-btn.external-btn:hover {
  background: #c0392b;
}

.update-tip {
  margin-top: 32px;
  padding: 16px;
  background: #f8f9fa;
  border-radius: 8px;
  border-left: 4px solid #3498db;
  text-align: center;
}

.update-tip p {
  margin: 0 0 8px 0;
  color: #555;
  font-size: 0.95rem;
}

.update-link {
  color: #3498db;
  text-decoration: none;
  font-weight: 500;
  transition: color 0.2s ease;
}

.update-link:hover {
  color: #2980b9;
  text-decoration: underline;
}

/* 响应式设计 */
@media (max-width: 768px) {
  .article-tree {
    padding: 16px;
  }

  .category-title {
    font-size: 1.5rem;
  }

  .article-header {
    flex-direction: column;
    align-items: flex-start;
    gap: 8px;
  }

  .article-footer {
    flex-direction: column;
    align-items: flex-start;
    gap: 12px;
  }

  .view-btn {
    width: 100%;
  }
}
</style>
